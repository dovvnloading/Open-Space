"""Check catalog triggers with two fixed sky apertures on original calibrated mosaics."""
import argparse
from http.client import IncompleteRead
import hashlib
import json
import math
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

from scipy import ndimage

from campaign import HERE, product_url

STEP = 0.04  # matched sky grid, arcsec
HALF = 51
YY, XX = np.mgrid[-HALF:HALF+1, -HALF:HALF+1]*STEP
SR_PER_ARCSEC2 = (np.pi/(180*3600))**2


def byte_range(url, start, size):
    for attempt in range(5):
        try:
            req = Request(url,headers={'Range':f'bytes={start}-{start+size-1}'})
            with urlopen(req,timeout=65) as response:
                expected=f'bytes {start}-{start+size-1}/'
                if response.status != 206 or not response.headers['Content-Range'].startswith(expected):
                    raise OSError('unexpected archive range')
                data=response.read()
                etag=response.headers['ETag'].strip('"')
            if len(data)!=size:
                raise OSError('truncated range')
            return data,dict(start=start,bytes=size,sha256=hashlib.sha256(data).hexdigest(),etag=etag)
        except (TimeoutError,OSError,IncompleteRead):
            if attempt==4:
                raise
            time.sleep(1+attempt)


def headers(product):
    url=product_url(product,'i2d.fits')
    data,receipt=byte_range(url,0,86400)
    pos=0
    while pos<len(data):
        start=pos
        while pos+80<=len(data) and data[pos:pos+8]!=b'END     ':
            pos+=80
        if pos+80>len(data):
            raise ValueError(f'incomplete FITS header for {product}')
        pos=2880*((pos+80+2879)//2880)
        h=fits.Header.fromstring(data[start:pos].decode('ascii'),sep='')
        if h.get('EXTNAME')=='SCI':
            if h['BUNIT']!='MJy/sr' or h['BITPIX']!=-32:
                raise ValueError(f'unexpected science plane in {product}')
            return url,h,pos,receipt
        if h['NAXIS']!=0:
            raise ValueError(f'unexpected HDU preceding science in {product}')
    raise ValueError(product)


def stamp(meta,ra,dec):
    url,h,offset,head_receipt=meta
    w=WCS(h,naxis=2)
    sx,sy=w.world_to_pixel_values(ra+XX/(3600*np.cos(np.deg2rad(dec))),dec+YY/3600)
    nx,ny=h['NAXIS1'],h['NAXIS2']
    y0=max(0,math.floor(sy.min())-1)
    y1=min(ny,math.ceil(sy.max())+2)
    if y1<=y0:
        raise ValueError('no stamp coverage')
    start=offset+4*y0*nx
    length=4*(y1-y0)*nx
    parts=[];receipts=[]
    for i in range(0,length,1024**2):
        data,r=byte_range(url,start+i,min(1024**2,length-i))
        if r['etag']!=head_receipt['etag']:
            raise ValueError('product version changed during download')
        parts.append(data);receipts.append(r)
    region=np.frombuffer(b''.join(parts),dtype='>f4').reshape(y1-y0,nx)
    pixels=ndimage.map_coordinates(region,[sy-y0,sx],order=1,mode='constant',cval=np.nan)
    return pixels,receipts


def aperture(image,x,y,radius):
    r=np.hypot(XX-x,YY-y)
    weights=np.clip((radius-r)/STEP+.5,0,1)
    sky=(r>radius+.12)&(r<radius+.32)
    if not np.isfinite(image[weights>0]).all() or np.isfinite(image[sky]).sum()<.9*sky.sum():
        return np.nan
    return np.sum((image[weights>0]-np.nanmedian(image[sky]))*weights[weights>0])*STEP**2*SR_PER_ARCSEC2*1e15


def measure(task,meta,cache):
    cid,band,ref,late,ra,dec,dx,dy,scale=task
    images=[];receipts=[]
    for p,r,d in [(ref,ra-dx/3600/np.cos(np.deg2rad(dec)),dec-dy/3600),(late,ra,dec)]:
        savename=cache/f'{cid}_{p}.npz'
        if savename.exists():
            with np.load(savename) as cached:
                pix=cached['pixels']; rec=json.loads(str(cached['receipts']))
        else:
            pix,rec=stamp(meta[p],r,d)
            np.savez_compressed(savename,pixels=pix,receipts=json.dumps(rec))
        images.append(pix)
        receipts.extend(dict(candidate=cid,product=p,**v) for v in rec)
    result=[]
    for radius in (.16,.28):
        a,b=[aperture(im,0,0,radius)*s for im,s in zip(images,(scale,1))]
        # Same apertures on the observed surroundings retain drizzle covariance,
        # local surface structure and inter-epoch subtraction residuals.
        background=[]
        for x in np.arange(-1.38,1.39,.46):
            for y in np.arange(-1.38,1.39,.46):
                if np.hypot(x,y)<.7:
                    continue
                background.append(aperture(images[1],x,y,radius)-scale*aperture(images[0],x,y,radius))
        bg=np.array(background)
        valid=np.isfinite(bg)
        noise=1.4826*np.nanmedian(abs(bg-np.nanmedian(bg))) if valid.sum()>=16 else np.nan
        err=max(noise,.025*(abs(a)+abs(b))/2) if np.isfinite(noise) else np.nan
        result.append(dict(candidate=cid,band=band,radius_arcsec=radius,ref_aperture_njy=a,search_aperture_njy=b,
                           delta_njy=b-a,empirical_err_njy=err,snr=(b-a)/err if err>0 else np.nan,blank_centers=int(valid.sum())))
    return result,receipts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache',type=Path)
    args=parser.parse_args();args.cache.mkdir(parents=True,exist_ok=True)
    triggers=pd.read_csv(HERE/'coincidences.csv')
    pairs=pd.read_csv(HERE/'pairs.csv').merge(pd.read_csv(HERE/'screen_counts.csv'),on=['pair','band'])
    tasks=[]
    for hit in triggers.itertuples():
        for band in hit.bands.split(','):
            pair=pairs[(pairs.pair==hit.pair)&(pairs.band==band)].iloc[0]
            tasks.append((hit.candidate,band,pair.reference,pair.search,hit.ra,hit.dec,pair.dx_arcsec,pair.dy_arcsec,pair.flux_scale))
    products=sorted({p for t in tasks for p in t[2:4]})
    with ThreadPoolExecutor(8) as pool:
        meta=dict(zip(products,pool.map(headers,products)))
        out=list(pool.map(lambda t:measure(t,meta,args.cache),tasks))
    pd.DataFrame([r for part,_ in out for r in part]).to_csv(HERE/'pixel_apertures.csv',index=False,float_format='%.9g')
    pd.DataFrame([dict(candidate='HEADER',product=p,**m[3]) for p,m in meta.items()]+[r for _,part in out for r in part]).to_csv(HERE/'pixel_inputs.csv',index=False)
    print(len(tasks),'position/filter checks completed')


if __name__=='__main__':
    main()
