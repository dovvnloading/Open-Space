"""Forced VIS/H ridge profiles at catalogue-selected positions, with local parallel controls."""
import argparse
import csv
import hashlib
import io
import json
import math
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates
from astropy.io import fits
from astropy.wcs import WCS

HERE=Path(__file__).resolve().parent
IBE='https://irsa.ipac.caltech.edu/ibe/data/euclid/q1/MER/'

def fetch(key,r,size):
    url=IBE+key+f"?center={float(r['ra']):.8f},{float(r['dec']):.8f}deg&size={size}arcsec&gzip=false"
    for i in range(3):
        try:
            with urlopen(url,timeout=70) as res: b=res.read()
            if not b.startswith(b'SIMPLE'):raise IOError('non-FITS response')
            break
        except OSError:
            if i==2:raise
            time.sleep(5*2**i)
    with fits.open(io.BytesIO(b)) as h:
        image=h[0].data.astype('f4'); w=WCS(h[0].header)
        x,y=w.world_to_pixel_values(float(r['ra']),float(r['dec']))
    return image,float(x),float(y),dict(url=url,sha256=hashlib.sha256(b).hexdigest(),bytes=len(b))

def measure(image,x,y,pa,length,width):
    # Same measurement on both bands at the VIS-selected axis; width is an observed scale.
    a=np.deg2rad(pa); t=np.linspace(-length/2,length/2,max(24,round(length)))
    u=np.array([-5.,-4.,-3.,-.45,0,.45,3.,4.,5.])*width
    valid=np.isfinite(image)&(image!=0)
    image=np.where(valid,image,0)
    def sample(shift):
        xx=x+np.cos(a)*t-np.sin(a)*(shift+u[:,None]);yy=y+np.sin(a)*t+np.cos(a)*(shift+u[:,None])
        ok=map_coordinates(valid.astype('f4'),[yy,xx],order=1,mode='constant')>.99
        pix=map_coordinates(image,[yy,xx],order=1,mode='constant')
        pix[~ok]=np.nan
        if np.mean(ok)<.94:return None
        lo,ctr,hi=[z.mean(axis=0) for z in np.split(pix,3,axis=0)]
        with np.errstate(invalid='ignore'):
            bins=np.array([[np.nanmean(v) for v in np.array_split(z,8)] for z in (ctr,lo,hi)])
        if not np.isfinite(bins).all():return None
        f=bins[0]-.5*(bins[1]+bins[2]);return float(f.mean()),int(np.sum(bins[0]>np.maximum(bins[1],bins[2]))),float(max(f.max(),0)/max(1e-20,np.maximum(f,0).sum()))
    value=sample(0)
    if value is None:return dict(measured=False)
    flank=[z[0] for shift in np.r_[-60:-12:3,15:61:3] if (z:=sample(shift))]
    if len(flank)<10:return dict(measured=False)
    med=np.median(flank);mad=1.4826*np.median(abs(np.array(flank)-med))
    if mad<=0:return dict(measured=False)
    return dict(measured=True,mean=float(value[0]),n_bins=value[1],brightest_bin_fraction=value[2],
                sky_positions=len(flank),offset_median=float(med),offset_mad=float(mad),contrast_mad=float((value[0]-med)/mad))

def thumb(image):
    m=np.isfinite(image)&(image!=0)
    if not m.any():return np.zeros((180,180),'u1')
    v=image[m]; q50,q84=np.percentile(v,[50,84]);sc=max(q84-q50,1e-10)
    arr=np.clip((np.arcsinh((image-q50)/(2*sc))+1.0)/4.0,0,1)
    arr=np.nan_to_num(arr)
    return np.array(Image.fromarray((255*arr[::-1]).astype('u1')).resize((180,180),Image.Resampling.LANCZOS))

def target(args):
    r,keys,cache=args;oid=str(r['object_id']);dest=cache/(oid+'.json')
    if dest.exists() and 'VIS_ccwx' in json.loads(dest.read_text()):return
    size=42
    try:
        row={k:r[k] for k in ['object_id','tileid','ra','dec','mer_position_angle','mer_det_quality_flag','mer_semimajor_axis','mer_ellipticity']}
        pa=float(r['mer_position_angle']);a=float(r['mer_semimajor_axis'])
        length=float(np.clip(np.sqrt(12)*a,35,200))
        row['forced_length_arcsec']=length*.1
        thumbs=[]
        for b,width in [('VIS',2.),('H',3.6)]:
            im,x,y,rec=fetch(keys[(str(r['tileid']),b)],r,size)
            row[b]=measure(im,x,y,pa+90,length,width)
            row[b+'_ccwx']=measure(im,x,y,pa,length,width)
            row[b+'_input']=rec
            thumbs.append(thumb(im))
        Image.fromarray(np.concatenate(thumbs,axis=1)).save(cache/(oid+'.jpg'),quality=85)
        dest.write_text(json.dumps(row,allow_nan=False)+'\n')
    except Exception as exc:
        dest.write_text(json.dumps(dict(object_id=oid,error=repr(exc)))+'\n')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--products',type=Path,default=HERE.parent/'euclid/q1_products.csv');p.add_argument('--tiles',type=Path,default=HERE.parent/'euclid/q1_tiles.csv')
    p.add_argument('--selected',type=Path,default=HERE/'hats_selected.csv')
    p.add_argument('--residual',action='store_true',help='also screen flagged deep-field rows whose axes miss the bright-star test')
    p.add_argument('--workers',type=int,default=8)
    a=p.parse_args();a.cache.mkdir(exist_ok=True,parents=True)
    tiles={r['tile'] for r in csv.DictReader(a.tiles.open()) if r['field']!='LDN1641'}
    keys={(r['tile'],r['band_or_product']):r['relative_to_mer'] for r in csv.DictReader(a.products.open())}
    align={}
    if a.residual:
        align={r['object_id']:r['aligned']=='True' for r in csv.DictReader((HERE/'gaia_alignment.csv').open())}
    rows=[r for r in csv.DictReader(a.selected.open()) if str(r['tileid']) in tiles and (
        int(r['mer_det_quality_flag'] or -1)==0 or (a.residual and not align[r['object_id']]))]
    with ThreadPoolExecutor(a.workers) as pool:list(pool.map(target,[(r,keys,a.cache) for r in rows]))
    print('targeted extragalactic objects',len(rows),flush=True)
    out=[json.loads((a.cache/(str(r['object_id'])+'.json')).read_text()) for r in rows]
    outname='targeted_with_residuals.jsonl' if a.residual else 'targeted_measurements.jsonl'
    (HERE/outname).write_text(''.join(json.dumps(v,separators=(',',':'))+'\n' for v in out))
if __name__=='__main__':main()
