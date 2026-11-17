"""Tie the native Spitzer mosaic WCS to unsaturated Gaia DR3 point sources."""
from pathlib import Path
import csv,json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from scipy.optimize import least_squares
from measure_l788 import offset,sample,HERE

def ties(image,wcs,epoch,template,center):
    records=[]
    with (HERE/'gaia-ties.csv').open() as f:
        for obj in csv.DictReader(f):
            ra,dec=offset(float(obj['ra']),float(obj['dec']),float(obj['pmra'])*(epoch-2016)/1000,float(obj['pmdec'])*(epoch-2016)/1000)
            x,y=np.array(wcs.world_to_pixel_values(ra,dec));h,w=image.shape
            if not (12<x<w-12 and 12<y<h-12) or np.hypot(x-center[0],y-center[1])<35:continue
            iy,ix=np.mgrid[round(y)-5:round(y)+6,round(x)-5:round(x)+6];v=image[iy,ix];ok=np.isfinite(v)&(np.hypot(ix-x,iy-y)<5)
            if ok.sum()<65:continue
            xx=ix[ok]-x;yy=iy[ok]-y;v=v[ok]; bg=np.median(v);v=v-bg
            def model(t):return t[2]*sample(template,18+xx-t[0],18+yy-t[1])+t[3]+t[4]*xx+t[5]*yy
            sol=least_squares(lambda t:model(t)-v,[0,0,max(v)*10,0,0,0],bounds=([-2,-2,0,-np.inf,-np.inf,-np.inf],[2,2,np.inf,np.inf,np.inf,np.inf]))
            dx,dy,amp=sol.x[:3];scatter=1.4826*np.median(np.abs(sol.fun-np.median(sol.fun)))
            if scatter and amp*template[18,18]/scatter>12 and np.hypot(dx,dy)<1.5:
                records.append([obj['source_id'],round(dx,3),round(dy,3),round(amp*template[18,18]/scatter,1)])
    if len(records)<5:raise ValueError('insufficient WCS ties')
    xy=np.array([r[1:3] for r in records]); med=np.median(xy,axis=0);dev=1.4826*np.median(abs(xy-med),axis=0)
    for r in records:r.append(bool(np.all(abs(r[1:3]-med)<np.maximum(3*dev,.2))))
    valid=np.array([r[1:3] for r in records if r[-1]])
    return np.median(valid,axis=0),np.std(valid,axis=0),records
