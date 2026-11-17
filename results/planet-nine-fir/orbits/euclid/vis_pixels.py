"""Extract real DR1 SCI/RMS/FLG data at conditional track positions from long VIS exposures."""
import csv
import hashlib
import io
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from time import sleep
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import numpy as np
from scipy.ndimage import map_coordinates

from intersections import HERE, p, spacecraft

BUCKET = 'https://nasa-irsa-euclid-dr1.s3.us-east-1.amazonaws.com/'
CACHE = HERE/'cache'
CACHE.mkdir(exist_ok=True)
REC = HERE/'receipts'
REC.mkdir(exist_ok=True)
# All calibrated VIS full-frame products used here have a 6-block primary HDU
# followed by 144 equal-sized groups. Validate headers/IDs, types, dimensions,
# offsets and byte ranges for each data window before using any pixel.
# McCracken et al. Q1/DR1 VIS flag guidance: object/stellar bits carry
# astrophysical signal. Reject the detector-invalid bits, preserving real objects.
BAD_MASK=sum(1 << b for b in (0,1,2,3,4,5,7,8,9,12,15,17,21,22))
BLOCK = 2880
QUAD = 50803200
HDU = QUAD//3
NX = 2048
NY = 2066
STAMP = 15
PIXSCALE = 0.1  # arcsec, common tangent grid
MODEL = p.family(json.loads((HERE.parent/'orbit_input.json').read_text()))
SC = spacecraft()
TREF = Time('2024-12-27T07:30:00', scale='utc')
RA0, DEC0 = p.sky(MODEL,TREF.tt.jd-2451545,SC(TREF.tdb.jd-2451545))


def read(url, label, region=None):
    hdr = {'Range':f'bytes={region[0]}-{region[1]}'} if region else {}
    dest = CACHE/(label+'.bin')
    if dest.exists():
        return dest.read_bytes()
    for attempt in range(3):
        try:
            with urlopen(Request(url,headers=hdr),timeout=90) as resp:
                cr = resp.headers.get('Content-Range')
                if region and (resp.status!=206 or not cr or not cr.startswith(f'bytes {region[0]}-{region[1]}/')):
                    raise ValueError(f'Unexpected HTTP byte range {cr} for {label}')
                data=resp.read()
                if region and len(data)!=region[1]-region[0]+1:
                    raise ValueError(f'Truncated {label}')
                dest.write_bytes(data)
                (REC/(label+'.json')).write_text(json.dumps(dict(
                    url=url, range=region, content_range=cr, bytes=len(data),
                    sha256=hashlib.sha256(data).hexdigest(), utc=datetime.now(timezone.utc).isoformat()),
                    indent=2)+'\n')
                return data
        except Exception:
            if attempt==2: raise
            sleep(1+attempt)


def qheader(url, key, start, name):
    data=read(url,key+'_hdr',(start,start+4*BLOCK-1))
    f=io.BytesIO(data)
    h=fits.Header.fromfile(f)
    assert h['EXTNAME']==name, (name,h['EXTNAME'])
    assert h['NAXIS1']==NX and h['NAXIS2']==NY
    return h, f.tell()


def jobs():
    rows=[r for r in csv.DictReader((HERE/'coverage.csv').open()) if int(r['tracks_in_polygon'])
          and 'VIS_SWL-DET-' in r['product'] and '-1-0000000_' in r['product']]
    out=[]
    for obs in sorted({r['obs_id'] for r in rows}):
        query=BUCKET+'?'+urlencode({'list-type':2,'prefix':f'dr1/VIS_QUAD/{obs}/'})
        xml=ET.fromstring(read(query,'vis_listing_'+obs))
        keys=[x.text for x in xml.iter() if x.tag.endswith('}Key')]
        for row in [r for r in rows if r['obs_id']==obs]:
            prod=row['product']
            exp=prod.split('DET-')[1].split('-0000000')[0]
            key=next(k for k in keys if k.endswith('.json') and exp in k)
            qc=json.loads(read(BUCKET+key,exp+'_qc'))
            for k,ra in qc.items():
                if not k.endswith('.QUAD_RA'): continue
                dec=qc[k[:-2]+'DEC']
                ds=np.hypot((RA0-ra)*np.cos(np.radians(dec)), DEC0-dec)*3600
                selected=np.flatnonzero(ds<160)
                if len(selected):
                    out.append((row,exp,k.removesuffix('.QUAD_RA'),selected))
    return out


def cut(job):
    row, exp, quad, ids=job
    cached=CACHE/(exp+'_'+quad+'_stamps.npz')
    if cached.exists():
        with np.load(cached) as saved:
            if 'bad_mask' in saved and int(saved['bad_mask'])==BAD_MASK: return str(cached)
    a,b=map(int,quad.split('.')[0].split('-'))
    idx=((a-1)*6+b-1)*4 + 'EFGH'.index(quad[-1])
    start=6*BLOCK + idx*QUAD
    url=BUCKET+f"dr1/VIS_QUAD/{row['obs_id']}/{row['product']}"
    stem=f'{exp}_{quad}'
    h, hbytes=qheader(url, stem+'SCI',start,quad+'.SCI')
    wcs=WCS(h)
    x,y=wcs.all_world2pix(RA0[ids],DEC0[ids],0)
    loc=(x>20)&(y>20)&(x<NX-21)&(y<NY-21)
    ids=ids[loc]
    if len(ids)==0: return None
    dd=np.arange(STAMP)-(STAMP-1)/2
    dx,dy=np.meshgrid(dd,dd)
    rr=RA0[ids,None,None]+dx*PIXSCALE/(3600*np.cos(np.radians(DEC0[ids,None,None])))
    dd=DEC0[ids,None,None]+dy*PIXSCALE/3600
    x,y=wcs.all_world2pix(rr,dd,0)
    lo=int(np.floor(y.min()))-2;hi=int(np.ceil(y.max()))+3
    assert 0<=lo<hi<=NY
    planes=[]
    for j,(name,dtype) in enumerate((('SCI','>f4'),('RMS','>f4'),('FLG','>i4'))):
        base=start+j*HDU
        if j:
            plane_h,hbytes=qheader(url,stem+name,base,quad+'.'+name)
            assert plane_h['BITPIX']==(-32 if name=='RMS' else 32)
        raw=read(url,stem+name+'_rows',(base+hbytes+lo*NX*4,base+hbytes+hi*NX*4-1))
        plane=np.frombuffer(raw,dtype=dtype).reshape(hi-lo,NX)
        # Require all four interpolation neighbours to have no instrument-invalid flags.
        if name=='FLG':
            plane=((plane & BAD_MASK)!=0).astype(float)
        planes.append(map_coordinates(plane.astype(float),[y-lo,x],order=1,mode='constant',cval=np.nan,prefilter=False).astype('f4'))
    np.savez(cached, bad_mask=BAD_MASK, ids=ids, sci=planes[0], rms=planes[1], flag=planes[2],
             utc=row['utc_mid'], exptime=float(h['EXPTIME']), zp=float(h['MAGZEROP']),
             extension=quad, product=row['product'])
    print(exp,quad,len(ids),flush=True)
    return str(cached)


if __name__=='__main__':
    tasks=jobs()
    print(f'{len(tasks)} candidate quadrants',flush=True)
    with ThreadPoolExecutor(max_workers=10) as pool:
        futs=[pool.submit(cut,j) for j in tasks]
        for f in as_completed(futs): f.result()
