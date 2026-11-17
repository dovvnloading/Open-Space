"""Resample downloaded flag-masked, calibrated VIS strips to a common stationary sky grid."""
import csv
import io
import json
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import map_coordinates
from scipy.spatial import ConvexHull

from vis_pixels import BAD_MASK, CACHE, HERE, MODEL, RA0, DEC0, SC, TREF, NX, NY, BLOCK, QUAD, HDU, PIXSCALE
from intersections import p

MEAN = np.array([RA0.mean(),DEC0.mean()])
COS = np.cos(np.radians(MEAN[1]))
STEP=PIXSCALE/3600


def xy(ra,de):
    return (ra-MEAN[0])*COS/STEP, (de-MEAN[1])/STEP


def world(x,y):
    return MEAN[0]+x*STEP/COS, MEAN[1]+y*STEP


def main():
    x,y=xy(RA0,DEC0)
    # Margin contains neighboring pixels, local PSFs and the model's motion.
    x0,x1=int(np.floor(x.min()))-20,int(np.ceil(x.max()))+21
    y0,y1=int(np.floor(y.min()))-20,int(np.ceil(y.max()))+21
    size=(y1-y0,x1-x0)
    out=CACHE/'mosaic'
    out.mkdir(exist_ok=True)
    rows=[]
    for old in out.glob("*.npy"): old.unlink()
    for fp in sorted(CACHE.glob('*_stamps.npz')):
        d=np.load(fp)
        quad=str(d['extension'])
        stem=fp.name.removesuffix('_stamps.npz')
        exp=stem.split('_')[0]
        a,b=map(int,quad.split('.')[0].split('-'))
        idx=((a-1)*6+b-1)*4+'EFGH'.index(quad[-1])
        # Use the validated and hashed FITS headers on the downloaded row strips.
        headers=[]; arrays=[]; starts=[]
        for k,typ in enumerate(['SCI','RMS','FLG']):
            base=6*BLOCK+QUAD*idx+k*HDU
            hbuf=io.BytesIO((CACHE/(stem+typ+'_hdr.bin')).read_bytes());h=fits.Header.fromfile(hbuf)
            receipt=json.loads((HERE/'receipts'/(stem+typ+'_rows.json')).read_text())
            start=(receipt['range'][0]-base-hbuf.tell())//(NX*4)
            raw=np.memmap(CACHE/(stem+typ+'_rows.bin'),dtype='>i4' if k==2 else '>f4',mode='r')
            arrays.append(raw.reshape(-1,NX));starts.append(start);headers.append(h)
        assert len(set(starts))==1
        lo=starts[0];hi=lo+arrays[0].shape[0]
        wcs=WCS(headers[0])
        nra,nde=wcs.all_pix2world([0,NX-1,NX-1,0],[lo,lo,hi-1,hi-1],0)
        xx,yy=xy(nra,nde)
        xmin=max(x0,int(np.floor(xx.min()))-2);xmax=min(x1,int(np.ceil(xx.max()))+3)
        ymin=max(y0,int(np.floor(yy.min()))-2);ymax=min(y1,int(np.ceil(yy.max()))+3)
        names=[out/(exp+'_'+s+'.npy') for s in ('data','weight')]
        mode='r+' if names[0].exists() else 'w+'
        maps=[np.lib.format.open_memmap(k,mode=mode,dtype='f4',shape=size) for k in names]
        # Fit a smooth local native-pixel mapping every 80 output pixels; check
        # against exact TPV WCS on interlaced midpoints before resampling.
        gx=np.linspace(xmin,xmax,max(4,(xmax-xmin)//80+2))
        gy=np.linspace(ymin,ymax,max(4,(ymax-ymin)//80+2))
        xx,yy=np.meshgrid(gx,gy)
        ux,uy=wcs.all_world2pix(*world(xx,yy),0)
        surf=[RectBivariateSpline(gy,gx,a) for a in (ux,uy)]
        tx=(gx[:-1]+gx[1:])/2;ty=(gy[:-1]+gy[1:])/2
        true=wcs.all_world2pix(*world(*np.meshgrid(tx,ty)),0)
        err=max(np.max(abs(sp(ty,tx)-t)) for sp,t in zip(surf,true))
        if err>0.005: raise ValueError(f'WCS interpolation error {err:.4f} pixels')
        counts=float(d['exptime'])*10**(.4*(float(d['zp'])-25))
        # Median sky on a deterministic detector subsample, measured separately
        # for each calibrated quadrant, leaving static source structure intact.
        sub=arrays[0][::16,::16]; good=(arrays[2][::16,::16]&BAD_MASK)==0
        if not np.any(good): continue
        bg=np.nanmedian(sub[good])
        for ystart in range(ymin,ymax,128):
            yend=min(ymax,ystart+128)
            u,v=(sp(np.arange(ystart,yend),np.arange(xmin,xmax)) for sp in surf)
            inside=(u>=2)&(u<NX-3)&(v>=lo+2)&(v<hi-3)
            coords=[v-lo,u]
            fg=map_coordinates(((arrays[2]&BAD_MASK)!=0).astype('f4'),coords,order=1,mode='constant',cval=1)
            rms=map_coordinates(arrays[1],coords,order=1,mode='constant',cval=0)
            data=map_coordinates(arrays[0],coords,order=1,mode='constant',cval=np.nan)
            good=inside&(fg==0)&(rms>0)&np.isfinite(data)
            region=np.s_[ystart-y0:yend-y0,xmin-x0:xmax-x0]
            maps[0][region][good]=(data[good]-bg)/counts
            maps[1][region][good]=(counts/rms[good])**2
        rows.append([stem,lo,hi,xmin,xmax,ymin,ymax,err,float(bg),counts])
        print(stem,round(err,7),flush=True)
    for m in maps: m.flush()
    with (HERE/'mosaic_wcs.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['quadrant','row0','row1','x0','x1','y0','y1','max_wcs_error_pix','bg_ADU','ADU_at_VIS25'])
        w.writerows(rows)
    (out/'geometry.json').write_text(json.dumps(dict(x0=x0,y0=y0,width=size[1],height=size[0],ra0=MEAN[0],dec0=MEAN[1],pix_as=PIXSCALE)))


if __name__=='__main__':
    main()
