"""Recover 6.1-arcsec Euclid VIS science/flag cutouts via verified HTTP byte ranges.

Run at repository root after imaging.py and archives.py. Astropy supplies the
per-quadrant TPV WCS. Fit no host model. Return repeatable field peaks only.
"""
import csv
import json
import re
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy.ndimage import gaussian_filter, maximum_filter, median_filter

from crossmatch import coords
from fetch import OUT, SHA, write_csv

CACHE=Path('.survey-cache/euclid_scene')
NAMES={'KMT-2024-BLG-0456':('17:50:26.78','-28:44:27.82'),'KMT-2024-BLG-3429':('17:52:00.03','-29:55:52.61'),
       'MOA-2025-BLG-0211':('17:55:00.6803','-30:33:42.6693')}
warnings.simplefilter('ignore',FITSFixedWarning)


def part(url,start,length,rec):
    file=CACHE/(url.split('/')[-1]+f'.{start}.{length}')
    info=Path(str(file)+'.json')
    if file.exists() and info.exists():
        data=file.read_bytes()
        entry=json.loads(info.read_text())
        assert SHA(data)==entry['sha256']
    else:
        request=Request(url,headers={'Range':f'bytes={start}-{start+length-1}'})
        with urlopen(request,timeout=30) as resp:
            match=re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)',resp.getheader('Content-Range',''))
            if resp.status!=206 or not match or int(match[1])!=start or int(match[2])+1!=start+length:
                raise ValueError('server did not return the specified FITS byte range')
            data=resp.read(length+1)
            if len(data)!=length:raise ValueError('wrong byte count')
        entry=dict(url=url,start=start,bytes=length,remote_file_bytes=int(match[3]),sha256=SHA(data))
        file.write_bytes(data);info.write_text(json.dumps(entry))
    rec.append(entry)
    return data


def header(raw):
    end=next(i+80 for i in range(0,len(raw),80) if raw[i:i+8].strip()==b'END')
    return fits.Header.fromstring(raw[:end].decode('ascii'),sep=''),(end+2879)//2880*2880


def cutout(url,ra,dec,rec):
    _,offset=header(part(url,0,23040,rec))
    stride,rem=divmod(rec[-1]['remote_file_bytes']-offset,144)
    if rem:raise ValueError('unexpected layout for 144 VIS quadrants')
    memo={}
    def sci(q):
        if q not in memo:
            begin=offset+q*stride
            h,size=header(part(url,begin,8640,rec))
            if not h['EXTNAME'].endswith('.SCI') or h['BITPIX']!=-32:
                raise ValueError('science plane layout mismatch')
            memo[q]=(h,begin+size)
        return memo[q]
    def center(q):
        h=sci(q)[0]
        return np.array([h['CRVAL1'],h['CRVAL2']])
    orig=center(0)
    dx=(center(20)-orig)/5;dy=(center(120)-orig)/5
    col,row=np.linalg.solve(np.array([dx,dy]).T,np.array([ra,dec])-orig)
    # Approximate CCD order guides header access; each final membership uses its own WCS.
    cells=sorted(range(36),key=lambda k:(k%6+.23-col)**2+(k//6-.22-row)**2)
    for cell in cells[:4]:
        for q in range(4*cell,4*cell+4):
            h,start=sci(q)
            w=WCS(h)
            x,y=w.all_world2pix(ra,dec,0)
            nx,ny=h['NAXIS1'],h['NAXIS2']
            if not (30<=x<nx-30 and 30<=y<ny-30):continue
            left,down=int(round(float(x)))-30,int(round(float(y)))-30
            def patch(start,bits):
                width=abs(bits)//8
                n=60*nx+61
                raw=part(url,start+(down*nx+left)*width,n*width,rec)
                typ='>f4' if bits==-32 else f'>i{width}'
                data=np.frombuffer(raw,typ)
                return np.stack([data[i*nx:i*nx+61] for i in range(61)])
            scene=patch(start,-32).astype(float)
            offset2=start+(nx*ny*4+2879)//2880*2880
            flaghead,flagbytes=header(part(url,offset2,8640,rec))
            if flaghead['EXTNAME'].endswith('.RMS'):
                rms_size=nx*ny*abs(flaghead['BITPIX'])//8
                offset2+=flagbytes+(rms_size+2879)//2880*2880
                flaghead,flagbytes=header(part(url,offset2,8640,rec))
            if not flaghead['EXTNAME'].endswith('.FLG') or flaghead['NAXIS1']!=nx:
                raise ValueError('unrecognized VIS flag extension '+flaghead['EXTNAME'])
            flag=patch(offset2+flagbytes,flaghead['BITPIX']).astype(np.int64)
            flag=flag*int(flaghead.get('BSCALE',1))+int(flaghead.get('BZERO',0))
            return h,scene,flag,w,left,down
    raise ValueError('no containing science patch at safe distance from quadrant edges')


def measure(job):
    event,ra,dec,p=job
    url=p['access_url']
    rec=[];detections=[]
    summary=dict(event=event,mjd=p['start_mjd'],url=url)
    try:
        h,im,flag,w,left,down=cutout(url,ra,dec,rec)
        smooth=gaussian_filter(im,.65)
        resid=smooth-median_filter(smooth,11)
        # VIS flag bits 18 and 24 identify STARSIGNAL and OBJECTS; both mark
        # astronomical light. Mask every other set bit, including bad pixels.
        clean=((flag & ~((1<<18)|(1<<24)))==0)&np.isfinite(resid)
        med=np.nanmedian(resid[clean]);mad=1.4826*np.nanmedian(np.abs(resid[clean]-med))
        yg,xg=np.indices(im.shape)
        ys,xs=np.where((resid==maximum_filter(resid,5))&(resid>med+6*mad)&clean)
        for y,x in zip(ys,xs):
            if min(y,x,60-y,60-x)<6:continue
            r=np.hypot(yg-y,xg-x)
            ann=(r>=4)&(r<=6)&clean
            if ann.sum()<16 or np.any(~clean[r<=1.5]):continue
            sky=np.median(im[ann])
            ap=(r<=2.3)&clean;core=(r<=1.8)&clean
            weights=np.maximum(im[core]-sky,0.)
            if weights.sum()==0:continue
            xx=np.average(xg[core],weights=weights);yy=np.average(yg[core],weights=weights)
            a,d=w.all_pix2world(xx+left,yy+down,0)
            east=(float(a)-ra)*np.cos(np.deg2rad(dec))*3600
            north=(float(d)-dec)*3600
            radius=float(np.hypot(east,north));flux=float(np.sum(im[ap]-sky))
            if radius>2 or flux<=0:continue
            detections.append(dict(event=event,mjd=p['start_mjd'],east_arcsec=round(east,5),north_arcsec=round(north,5),
                        radius_arcsec=round(radius,5),net_aperture_counts=round(flux,3),filtered_peak_over_mad=round(float((resid[y,x]-med)/mad),2)))
        summary.update(status='ok',extension=h['EXTNAME'],patch_left_pixel=left,patch_bottom_pixel=down,usable_pixels=int(clean.sum()),total_pixels=im.size,
                       detected_within_two_arcsec=len(detections),local_filtered_mad=round(float(mad),3))
    except Exception as e:
        summary.update(status='error',error=repr(e))
    return summary,detections,[dict(event=event,**r) for r in rec]


def clusters(detections, names=None):
    out=[]
    for event in (NAMES if names is None else names):
        groups=[]
        for row in sorted((r for r in detections if r['event']==event),key=lambda r:(r['mjd'],r['radius_arcsec'])):
            point=np.array([row['east_arcsec'],row['north_arcsec']])
            choices=[(np.linalg.norm(point-np.median([m[0] for m in g],axis=0)),j) for j,g in enumerate(groups)
                     if row['mjd'] not in [m[1]['mjd'] for m in g]]
            if choices and min(choices)[0]<.10:groups[min(choices)[1]].append((point,row))
            else:groups.append([(point,row)])
        for g in groups:
            pos=np.array([m[0] for m in g]);flux=np.array([m[1]['net_aperture_counts'] for m in g]);loc=np.median(pos,axis=0)
            out.append(dict(event=event,independent_dithers=len(g),east_arcsec=round(loc[0],4),north_arcsec=round(loc[1],4),
                radius_arcsec=round(float(np.linalg.norm(loc)),4),rms_offset_arcsec=round(float(np.mean(np.sum((pos-loc)**2,axis=1))**.5),4),
                median_counts=round(float(np.median(flux)),2),flux_max_over_min=round(float(max(flux)/min(flux)),3)))
    return sorted(out,key=lambda r:(r['event'],r['radius_arcsec']))


def repeated_pairs(detections,sources,names=None):
    rows=[]
    for event in (NAMES if names is None else names):
        repeats=[s for s in sources if s['event']==event and s['independent_dithers']>=3][:2]
        data=[]
        dates=sorted({r['mjd'] for r in detections if r['event']==event})
        for epoch in dates:
            matches=[]
            for s in repeats:
                found=[r for r in detections if r['event']==event and r['mjd']==epoch
                       and np.hypot(r['east_arcsec']-s['east_arcsec'],r['north_arcsec']-s['north_arcsec'])<.1]
                if found:matches.append(found[0])
            if len(matches)!=2:continue
            a,b=matches
            sep=np.hypot(a['east_arcsec']-b['east_arcsec'],a['north_arcsec']-b['north_arcsec'])
            dm=2.5*np.log10(a['net_aperture_counts']/b['net_aperture_counts'])
            data.append((sep,dm))
        if data:
            (lo,dl),(mid,dm),(hi,dh)=np.quantile(data,[0,.5,1],axis=0)
            rows.append(dict(event=event,paired_exposures=len(data),median_pair_sep_arcsec=round(mid,4),
                             min_sep_arcsec=round(lo,4),max_sep_arcsec=round(hi,4),
                             approximate_delta_vis_mag=round(dm,3),min_delta_mag=round(dl,3),max_delta_mag=round(dh,3)))
    return rows


if __name__=='__main__':
    CACHE.mkdir(parents=True,exist_ok=True)
    inv=list(csv.DictReader((OUT/'imaging_products.csv').open()))
    for block in json.loads((OUT.parent/'archive-queries.json').read_text()):inv.extend(block['products'])
    jobs=[]
    for name,sex in NAMES.items():
        ra,dec=coords(*sex)
        urls={p['access_url']:p for p in inv if p['event']==name and '-DET-' in (p.get('access_url') or '')}
        products=sorted(urls.values(),key=lambda p:float(p['start_mjd']))
        assert len(products)==16,(name,len(products))
        jobs.extend((name,ra,dec,products[i]) for i in (0,1,2,3,5,9,13))
    with ThreadPoolExecutor(8) as pool:
        results=list(pool.map(measure,jobs))
    det=[r for _,rows,_ in results for r in rows]
    write_csv(OUT/'euclid_scene_checks.csv',[s for s,_,_ in results])
    write_csv(OUT/'euclid_scene_receipts.csv',[r for _,_,rows in results for r in rows])
    write_csv(OUT/'euclid_scene_detections.csv',det)
    sources=clusters(det)
    write_csv(OUT/'euclid_scene_sources.csv',sources)
    write_csv(OUT/'euclid_scene_pairs.csv',repeated_pairs(det,sources))
