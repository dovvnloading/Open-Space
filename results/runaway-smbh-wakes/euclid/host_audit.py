"""Compare Euclid VIS ridge endpoints with the Q1 MER catalogue on the sky.

Coordinates in detector JSON refer to zero-origin native pixels. Positives and
sign-inverted controls use exactly the same angular host rule. A match is a
geometric association only. Raw TAP replies can be regenerated from the URLs.
"""
import argparse
import csv
import hashlib
import gzip
import io
import json
import math
import time
from xml.etree import ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from acquire import write_csv

ESA = 'https://eas.esac.esa.int/tap-server/tap/sync'
MER_COLS = ('object_id right_ascension declination extended_flag extended_prob point_like_prob '
            'gaia_id vis_det semimajor_axis ellipticity position_angle flux_vis_2fwhm_aper det_quality_flag').split()


def tap(endpoint, query, directory):
    url = endpoint + '?' + urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', MAXREC=10000, QUERY=query))
    key = hashlib.sha256(url.encode()).hexdigest()
    dest = directory / (key+'.csv')
    if not dest.exists():
        for attempt in range(3):
            try:
                with urlopen(url, timeout=90) as response:
                    data = response.read()
                if not data or data.lstrip()[:1] == b'<':
                    raise ValueError(data[:400].decode(errors='replace'))
                tmp = dest.with_suffix('.tmp')
                tmp.write_bytes(data)
                tmp.replace(dest)
                break
            except (OSError, ValueError):
                if attempt == 2:
                    raise
                time.sleep(2 ** attempt)
    data = dest.read_bytes()
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    if len(rows) >= 10000:
        raise ValueError('TAP response may have reached MAXREC')
    return rows, dict(query=query, url=url, sha256=hashlib.sha256(data).hexdigest(),
                      bytes=len(data), rows=len(rows))


def gator_mer(tile, ra, dec, radius_arcsec, directory):
    cols = ['ra' if c == 'right_ascension' else 'dec' if c == 'declination' else c for c in MER_COLS]
    url = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?' + urlencode(dict(
        catalog='euclid_q1_mer_catalogue', spatial='Cone', radius=f'{radius_arcsec:.6f}', radunits='arcsec',
        objstr=f'{ra:.9f} {dec:.9f}', outfmt=3, selcols=','.join(cols)))
    key = hashlib.sha256(url.encode()).hexdigest()
    path = directory / (key + '.xml')
    if not path.exists():
        for attempt in range(3):
            try:
                with urlopen(url, timeout=45) as res:
                    raw = res.read()
                tree = ET.fromstring(raw)
                if tree.find('.//{*}TABLEDATA') is None:
                    raise ValueError('Gator result has no TABLEDATA: '+raw[:140].decode(errors='replace'))
                path.write_bytes(raw)
                break
            except (ValueError,OSError,ET.ParseError) as exc:
                if attempt==2: raise
                delay=2**attempt
                if getattr(exc,'code',None)==429:
                    retry=exc.headers.get('Retry-After','60')
                    delay=int(retry) if retry.isdigit() else 60
                time.sleep(delay)
    raw=path.read_bytes();tree=ET.fromstring(raw)
    fs=[e.get('name') for e in tree.findall('.//{*}FIELD')]
    rows=[dict(zip(fs,[e.text or '' for e in tr.findall('{*}TD')])) for tr in tree.findall('.//{*}TR')]
    if not set(cols)<=set(fs):raise ValueError('Gator returned unexpected fields')
    for row in rows:
        row['right_ascension']=row['ra'];row['declination']=row['dec']
        for key in list(row):
            if row[key].lower() in ('null','nan'): row[key]=''
    return rows,dict(url=url,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),rows=len(rows))


def angular_distance(point, ra, dec):
    return 3600 * np.hypot((point[0]-ra)*np.cos(np.deg2rad(point[1])), point[1]-dec)


def number(value):
    return float(value) if value and value!='null' else np.nan


def match(tile, ridge, wcs, objects, sign, n, bright_stars=None):
    center, e0, e1 = wcs.all_pix2world([(ridge['x'],ridge['y']),
                                       (ridge['x0'],ridge['y0']), (ridge['x1'],ridge['y1'])], 0)
    ra = np.array([float(r['right_ascension']) for r in objects])
    dec = np.array([float(r['declination']) for r in objects])
    d = np.array([angular_distance(end,ra,dec) for end in (e0,e1)])
    # Q1's extended_flag/extended_prob are empty in these cones. Use the
    # populated VIS point-source probability and size as a host proxy.
    ext = np.array([r['vis_det']=='1' and not r['gaia_id'] and
                    number(r['point_like_prob'])<0.5 and number(r['semimajor_axis'])>=2.
                    for r in objects], dtype=bool)
    hosts = np.where(ext & (np.min(d,axis=0)<=1.0))[0]
    out = dict(id=f'{tile}_{sign}_{n}', tile=tile, sign=sign, index=n, x0=ridge['x0'], y0=ridge['y0'],
               x1=ridge['x1'], y1=ridge['y1'], ra_deg=center[0], dec_deg=center[1],
               ra0_deg=e0[0], dec0_deg=e0[1], ra1_deg=e1[0], dec1_deg=e1[1], length_arcsec=ridge['length_arcsec'],
               axis_ratio=ridge['axis_ratio'], score=ridge['score'], distinct_extended_hosts=int(len(hosts)),
               host_end='', host_object_id='', host_sep_arcsec='', host_ra='', host_dec='', host_vis_flux_uJy='',
               nearest_gaia_id='', nearest_gaia_arcsec='', nearest_gaia_vis_flux_uJy='', nearest_any_arcsec=float(d.min()) if len(objects) else '')
    if len(hosts):
        i=int(hosts[np.argmin(d[:,hosts].min(axis=0))]);end=int(d[:,i].argmin());h=objects[i]
        out.update(host_end=end,host_object_id=h['object_id'],host_sep_arcsec=d[end,i],host_ra=ra[i],host_dec=dec[i],
                   host_vis_flux_uJy=h['flux_vis_2fwhm_aper'])
    dx=(ra-center[0])*np.cos(np.deg2rad(center[1]))*3600
    dy=(dec-center[1])*3600
    delta=np.array([(e1[0]-e0[0])*np.cos(np.deg2rad(center[1])),e1[1]-e0[1]])*3600
    length=np.linalg.norm(delta);v=delta/length
    along=dx*v[0]+dy*v[1];cross=-dx*v[1]+dy*v[0]
    inside=(np.abs(along)<length/2)&(np.abs(cross)<=0.6)
    out['interior_mer_ids']=' '.join(objects[i]['object_id'] for i in np.flatnonzero(inside))
    out['interior_offsets_arcsec']=' '.join(f'{along[i]:.3f},{cross[i]:.3f}' for i in np.flatnonzero(inside))
    def axis_host(along,cross):
        distance=np.abs(along)-length/2
        eligible=ext & (distance>=0) & (distance<=5) & (np.abs(cross)<=.6)
        candidates=np.flatnonzero(eligible)
        if not len(candidates): return '', '', '', ''
        i=min(candidates,key=lambda k:np.hypot(distance[k],cross[k]))
        return objects[i]['object_id'],float(distance[i]),float(cross[i]),int(along[i]>0)
    out['axis_host_id'],out['axis_gap_arcsec'],out['axis_offset_arcsec'],out['axis_host_end']=axis_host(along,cross)
    axis_match = next((r for r in objects if r['object_id']==out['axis_host_id']), {})
    out['axis_host_ra']=axis_match.get('right_ascension','');out['axis_host_dec']=axis_match.get('declination','')
    out['rot90_axis_host_id'],_,_,_=axis_host(cross,-along)
    out['rot90_core_host']=int(np.any(ext & (np.hypot(np.abs(cross)-length/2,along)<=1.)))
    out.update(aligned_gaia_id='',aligned_gaia_axis_arcsec='',aligned_gaia_cross_arcsec='',aligned_gaia_vis_flux_uJy='')
    stars=[i for i,r in enumerate(objects) if r['gaia_id'] not in ('','0','-1','null','NaN')]
    if stars:
        ds=angular_distance(center,ra,dec);i=min(stars,key=lambda k:ds[k]);h=objects[i]
        out.update(nearest_gaia_id=h['gaia_id'], nearest_gaia_arcsec=ds[i],
                   nearest_gaia_vis_flux_uJy=h['flux_vis_2fwhm_aper'])
        aligned=[k for k in stars if abs(cross[k])<=1.2 and abs(along[k])>length/2 and number(objects[k]['flux_vis_2fwhm_aper'])>3.]
        if aligned:
            i=min(aligned,key=lambda k:abs(cross[k]))
            out.update(aligned_gaia_id=objects[i]['gaia_id'],aligned_gaia_axis_arcsec=along[i],
                       aligned_gaia_cross_arcsec=cross[i],aligned_gaia_vis_flux_uJy=objects[i]['flux_vis_2fwhm_aper'])
    out.update(gaia_dr3_spike_id='',gaia_dr3_g='',gaia_dr3_axis_arcsec='',gaia_dr3_cross_arcsec='')
    if bright_stars is not None:
        sx=(bright_stars['ra']-center[0])*np.cos(np.deg2rad(center[1]))*3600
        sy=(bright_stars['dec']-center[1])*3600
        sa=sx*v[0]+sy*v[1];sc=-sx*v[1]+sy*v[0];sr=np.hypot(sx,sy)
        indx=np.flatnonzero((sr>=3)&(sr<=120)&(np.abs(sc)<=1)&(bright_stars['g']<16))
        if len(indx):
            k=indx[np.argmin(np.abs(sc[indx]))]
            out.update(gaia_dr3_spike_id=str(bright_stars['id'][k]),gaia_dr3_g=bright_stars['g'][k],
                       gaia_dr3_axis_arcsec=sa[k],gaia_dr3_cross_arcsec=sc[k])
    return out


def process(tile, detections, args):
    with fits.open(args.images/(tile+'_VIS.fits'), memmap=False) as hdus:
        h=hdus[0].header
        wcs=WCS(h)
    xy=(np.array([h['NAXIS1'],h['NAXIS2']])-1)/2
    ra,dec=wcs.all_pix2world([xy],0)[0]
    rad=(np.hypot(*xy)*.1+16)/3600
    objects,receipt=gator_mer(tile,ra,dec,rad*3600,args.cache)
    stars=None
    if getattr(args,'gaia_fields',None):
        field='EDF-N' if dec>55 else 'EDF-S' if dec < -40 else 'EDF-F'
        stars=args.gaia_fields[field]
    rows=[match(tile,d,wcs,objects,s,n,stars) for s,n,d in detections]
    return rows, dict(tile=tile, **receipt)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--search', type=Path, required=True, nargs='+')
    p.add_argument('--images', type=Path, required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--jobs',default=4,type=int)
    p.add_argument('--gaia-dir',type=Path,help='Folder containing gaia_EDF-{N,S,F}.csv.gz field extracts for radial-spike diagnostics')
    p.add_argument('--prefetch-tiles',type=Path,help='Fetch all extragalactic tile-center catalogs before detector completion')
    args=p.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True);args.out.mkdir(parents=True,exist_ok=True)
    args.gaia_fields={}
    if args.gaia_dir:
        for field in ('EDF-N','EDF-S','EDF-F'):
            with gzip.open(args.gaia_dir/('gaia_'+field+'.csv.gz'),'rt') as f:rows=list(csv.DictReader(f))
            args.gaia_fields[field]={k:np.array([r[col] for r in rows],dtype=dt) for k,col,dt in
                [('id','source_id',str),('ra','ra',float),('dec','dec',float),('g','phot_g_mean_mag',float)]}
    groups={}
    for directory in args.search:
        for path in sorted(directory.glob('*.json')):
            if not path.stem.split('_')[0].isdigit(): continue
            a=json.loads(path.read_text())
            if 'detections' not in a: continue
            tile=a['inputs']['science']['name'].split('_')[0]
            sign='negative' if a['params']['invert'] else 'positive'
            group=groups.setdefault(tile,{})
            for i,d in enumerate(a['detections']):
                if d['profile_gate']: group[sign,i]=d
    if args.prefetch_tiles:
        for r in csv.DictReader(args.prefetch_tiles.open()):
            if r['field']!='LDN1641': groups.setdefault(r['tile'],{})
    def run(tile):
        try:
            return process(tile,[(s,i,d) for (s,i),d in groups[tile].items()],args)
        except (OSError,ValueError,ET.ParseError) as exc:
            return [],dict(tile=tile,error=repr(exc))
    tiles=sorted(k for k,v in groups.items() if v or args.prefetch_tiles)
    with ThreadPoolExecutor(args.jobs) as pool:
        bundles=list(pool.map(run,tiles))
    rows=[r for group,_ in bundles for r in group]
    if rows: write_csv(args.out/'endpoint_candidates.csv',rows)
    (args.out/'host_queries.json').write_text(json.dumps([b for _,b in bundles],indent=1)+'\n')


if __name__=='__main__':
    main()
