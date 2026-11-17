"""Fetch FITS bytes, verify object identity, and retain small WCS-centered cutouts.

Usage: python read_cubes.py CACHE [--positions positions.csv] [--workers 3]
Temporary stamps are written to CACHE. Receipts reference complete immutable-by-hash inputs.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]/'neowise-nuclear-flares'))
from spherex_tools import tile_url

RADIUS=17

def unpack(data, rows, cache, tile, sink=None):
    with fits.open(io.BytesIO(data),memmap=False,do_not_scale_image_data=True) as hdus:
        img=hdus['IMAGE']
        assert img.header['NPIX']==tile and img.header['ORDER']==5 and img.header['BUNIT']=='MJy/sr'
        assert img.data.shape==(102,512,512) and hdus['NHITS'].data.shape==hdus['FLAGS'].data.shape==img.data.shape
        # The spectral -TAB axis is tabulated separately. Construct a spatial-only WCS.
        hdr=img.header.copy()
        for key in list(hdr):
            if key in ('NAXIS3','CTYPE3','CNAME3','CUNIT3','CRPIX3','CDELT3','CRVAL3','PC3_3') or key.startswith(('PS3_','PV3_')):
                del hdr[key]
        hdr['NAXIS']=hdr['WCSAXES']=2
        w=WCS(hdr)
        w.wcs.bounds_check(False,False)  # valid tiles can cross the 180-degree HPX seam
        meta={'start':hdus[0].header['DATE-BEG'],'end':hdus[0].header['DATE-END'],
              'pipeline':hdus[0].header['VERSION'],'positions':[]}
        for row in rows:
            x,y=map(float,w.world_to_pixel_values(float(row['ra_2025_9']),float(row['dec_2025_9'])))
            if not(-.5<=x<511.5 and -.5<=y<511.5):
                # HPX's longitude discontinuity at RA=180 has two equivalent projected X values.
                wrap=np.linalg.solve(w.pixel_scale_matrix,[360.,0.])
                solutions=[(x+sign*wrap[0],y+sign*wrap[1]) for sign in (-1,1)]
                inside=[(xx,yy) for xx,yy in solutions if -.5<=xx<511.5 and -.5<=yy<511.5]
                if len(inside)!=1:
                    raise ValueError(f"position {row['label']} projects outside its HEALPix tile: {(x,y)}")
                x,y=map(float,inside[0])
            lon,lat=w.pixel_to_world_values(x,y)
            if not np.isfinite(lon+lat) or abs((lon-float(row['ra_2025_9'])+180)%360-180)>1e-6 or abs(lat-float(row['dec_2025_9']))>1e-6:
                raise ValueError('HPX position failed sky-coordinate round trip')
            left,right=max(0,round(x)-RADIUS),min(512,round(x)+RADIUS+1)
            lower,upper=max(0,round(y)-RADIUS),min(512,round(y)+RADIUS+1)
            stamp={}
            for name in ('IMAGE','NHITS','FLAGS'):
                hdu=hdus[name]
                raw=hdu.data[:,lower:upper,left:right]
                # Astropy raw data retains the FITS signed-storage offset.
                if 'BZERO' in hdu.header or 'BSCALE' in hdu.header:
                    raw=raw.astype(np.int64)*hdu.header.get('BSCALE',1)+hdu.header.get('BZERO',0)
                stamp[name]=raw.astype({'IMAGE':'f4','NHITS':'i4','FLAGS':'u4'}[name])
            for name in ('WAVELENGTH','WL_MIN','WL_MAX','DETECTOR','SUBCHAN'):
                stamp[name]=hdus['SPECTRAL_CHANNELS'].data[name]
            stamp.update(wcs=hdr.tostring(),xy=np.array([x,y]),origin=[left,lower])
            if sink is None:
                np.savez_compressed(cache/(row['label']+'.npz'),**stamp)
            else:
                sink(row,stamp)
            meta['positions'].append(dict(label=row['label'],pixel_x=x,pixel_y=y))
        return meta

def fetch_tile(task, sink=None):
    tile,rows,cache,out=task
    url=tile_url(tile)
    marker=out/f'{tile}.json'
    if sink is None and marker.exists() and all((cache/(r['label']+'.npz')).exists() for r in rows): return
    for attempt in range(4):
        try:
            # Several positions share each object. Read whole objects in auditable
            # sequential HTTP byte ranges (including header and ancillary HDUs).
            start=datetime.now(timezone.utc).isoformat()
            data=bytearray()
            etag=None
            total=None
            step=24*1024**2
            while total is None or len(data)<total:
                lo=len(data)
                hi=lo+step-1 if total is None else min(lo+step,total)-1
                headers={'Range':f'bytes={lo}-{hi}'}
                if etag: headers['If-Match']=etag
                with urlopen(Request(url,headers=headers),timeout=120) as r:
                    bounds,slength=r.headers['Content-Range'].split('/')
                    if r.status!=206 or bounds!=f'bytes {lo}-{hi}': raise IOError(f'range mismatch {bounds}')
                    total=int(slength)
                    if etag and r.headers['ETag']!=etag: raise IOError('object changed')
                    etag=r.headers['ETag']
                    chunk=r.read()
                    if len(chunk)!=hi-lo+1: raise IOError('truncated HTTP response')
                    data.extend(chunk)
            if '-' not in etag and hashlib.md5(data).hexdigest()!=etag.strip('"'):
                raise IOError('archive ETag differs from MD5')
            meta=unpack(data,rows,cache,tile,sink=sink)
            receipt=dict(tile=tile,url=url,etag=etag,size=total,sha256=hashlib.sha256(data).hexdigest(),
                         retrieved_utc=start,**meta)
            marker.write_text(json.dumps(receipt,indent=2)+'\n')
            print(tile,'ok',flush=True)
            return
        except Exception as exc:
            print(tile,repr(exc),flush=True)
            if attempt==3: raise
            time.sleep(2**attempt)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('cache',type=Path)
    p.add_argument('--positions',type=Path,default=ROOT/'positions.csv')
    p.add_argument('--workers',type=int,default=3)
    a=p.parse_args()
    a.cache.mkdir(parents=True,exist_ok=True)
    out=ROOT/'cube_receipts'
    out.mkdir(exist_ok=True)
    groups=defaultdict(list)
    for row in csv.DictReader(a.positions.open()): groups[int(row['tile'])].append(row)
    with ThreadPoolExecutor(a.workers) as pool: list(pool.map(fetch_tile,[(t,rows,a.cache,out) for t,rows in groups.items()]))
if __name__=='__main__': main()
