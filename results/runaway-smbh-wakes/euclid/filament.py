"""Measure independent longitudinal segments of the central 70% of a trial line.

Translated copies of the identical aperture measure the local, correlated
noise. NISP/ VIS images must carry a calibrated MAGZERO. Screen all geometric
host matches, keeping every row irrespective of the follow-up result.
"""
import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

from acquire import IBE, fetch, sha, write_csv


def measure(a, p0, p1, scale, mask=None):
    delta = p1-p0
    length = np.linalg.norm(delta)
    v = delta / length
    normal = np.array([-v[1],v[0]])
    # The five equal segments span 15%-85% of the ridge, excluding endpoint galaxies.
    t = np.linspace(.15,.85,50,endpoint=False)+.007
    xy = p0 + t[:,None]*delta
    transverse = np.array([-.10,0,.10,-.70,-.60,-.50,.50,.60,.70])/scale
    traces = xy[:,None,:] + transverse[None,:,None]*normal
    offsets = np.array([(x,y) for x in (-5.,-3.,-1.,1.,3.,5.) for y in (-5.,-3.,-1.,1.,3.,5.)
                        if np.hypot(x,y)>2. and abs(np.dot((x,y),normal))>1.3])/scale
    offsets = np.vstack([[0.,0.],offsets])
    pts = traces[None,:,:,:]+offsets[:,None,None,:]
    xx,yy=pts[...,0],pts[...,1]
    values=map_coordinates(a,[yy.ravel(),xx.ravel()],order=1,mode='constant',cval=np.nan).reshape(xx.shape)
    if mask is not None:
        valid=map_coordinates(mask.astype('u1'),[yy.ravel(),xx.ravel()],order=1,mode='constant',cval=0).reshape(xx.shape)>0
        values[~valid]=np.nan
    # Require centreline above both parallel sides. Median along each segment
    # prevents a single bright knot from supplying most of the line signal.
    ridge=values[:,:,:3].mean(axis=2)-np.maximum(values[:,:,3:6].mean(axis=2),values[:,:,6:].mean(axis=2))
    binned=np.median(ridge.reshape(len(offsets),5,10),axis=2)
    good=np.isfinite(binned).all(axis=1)
    if not good[0] or good.sum()<12:
        return dict(usable=0, blank_apertures=int(good[1:].sum()),median_sigma='',positive_bins='',bins_above_2sigma='',
                    bins_sigma='',central_ridge_contrast_sigma='')
    blanks=binned[1:][good[1:]]
    bg=np.median(blanks,axis=0)
    sc=1.4826*np.median(np.abs(blanks-bg),axis=0)
    v0=np.median(blanks,axis=1);s0=1.4826*np.median(np.abs(v0-np.median(v0)))
    if np.any(sc<=0) or s0<=0:
        return dict(usable=0,blank_apertures=len(blanks),median_sigma='',positive_bins='',bins_above_2sigma='',
                    bins_sigma='',central_ridge_contrast_sigma='')
    z=(binned[0]-bg)/sc
    return dict(usable=1,blank_apertures=len(blanks),median_sigma=float(np.median(z)),positive_bins=int(np.sum(z>0)),
                bins_above_2sigma=int(np.sum(z>=2)),bins_sigma=' '.join(f'{v:.3f}' for v in z),
                central_ridge_contrast_sigma=float((np.median(binned[0])-np.median(v0))/s0))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ('candidates','products','images','cache','out'):p.add_argument('--'+flag,type=Path,required=True)
    p.add_argument('--jobs',type=int,default=6);p.add_argument('--nir',action='store_true')
    args=p.parse_args();args.cache.mkdir(exist_ok=True,parents=True);args.out.mkdir(exist_ok=True,parents=True)
    hosts=[r for r in csv.DictReader(args.candidates.open()) if int(r['distinct_extended_hosts']) or r['axis_host_id']]
    products={(r['tile'],r['band_or_product']):r['relative_to_mer'] for r in csv.DictReader(args.products.open())}
    def run(task):
        row,band=task;tile=row['tile'];receipt=None
        path=args.images/(tile+'_'+band+'.fits')
        if band=='VIS':
            with fits.open(args.images/(tile+'_rms.fits')) as hdus: rms=hdus[0].data
            mask=np.isfinite(rms)&(rms>0)&(rms<1e8)
        else:
            mask=None
        if not path.exists():
            url=IBE+'q1/MER/'+products[tile,band]+f"?center={float(row['ra_deg']):.8f},{float(row['dec_deg']):.8f}deg&size=24arcsec&gzip=false"
            path=args.cache/(row['id']+'_'+band+'.fits')
            raw=fetch(url,path)
            receipt=dict(id=row['id'],tile=tile,band=band,url=url,bytes=len(raw),sha256=sha(raw))
        with fits.open(path) as hdus:
            h=hdus[0].header;a=hdus[0].data;w=WCS(h)
            p0,p1=w.all_world2pix([[float(row['ra0_deg']),float(row['dec0_deg'])],
                                    [float(row['ra1_deg']),float(row['dec1_deg'])]],0)
            mask=(np.isfinite(a)&(a!=0)) if mask is None else (mask & np.isfinite(a) & (a!=0))
            r=measure(-a if row['sign']=='negative' and band=='VIS' else a,p0,p1,np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600,mask)
        return dict(id=row['id'],tile=tile,sign=row['sign'],band=band,**r),receipt
    first=list(map(run,[(row,'VIS') for row in hosts]))
    metrics=[r for r,_ in first];receipts=[]
    if metrics:write_csv(args.out/'filament_photometry.csv',metrics)
    if args.nir:
        # The VIS coherence cut is identical on both signs. Only sky positives get NISP follow-up.
        ids={r['id'] for r in metrics if r['usable'] and r['positive_bins']>=4 and r['central_ridge_contrast_sigma']>=3.}
        tasks=[(r,b) for r in hosts if r['id'] in ids and r['sign']=='positive' for b in ('Y','J','H')]
        with ThreadPoolExecutor(args.jobs) as pool:
            def measured(task):
                try: return run(task)
                except (OSError,ValueError) as exc:
                    row,band=task
                    return None,dict(id=row['id'],tile=row['tile'],band=band,error=repr(exc))
            for metric,receipt in pool.map(measured,tasks):
                if metric:metrics.append(metric)
                if receipt: receipts.append(receipt)
    if metrics:write_csv(args.out/'filament_photometry.csv',metrics)
    (args.out/'nir_queries.json').write_text(json.dumps(receipts,indent=1)+'\n')


if __name__=='__main__': main()
