"""Inject empirical companions around the nearest repeated Euclid VIS peak.

Test a hypothetical source identification. Each trial must create a new
resolved peak at its sky position in at least half of the usable frames.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter, median_filter, shift

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'other'))
from euclid_scene import cutout, coords, CACHE

RADII=[.2,.4,.7,1.]
DELTAS=[0.,1.,2.,3.,4.]
NANGLE=24


def peaks(residual, valid, threshold):
    return np.transpose(np.nonzero(valid & (residual>threshold) & (maximum_filter(residual,5)==residual)))


def nearby(positions,x,y):
    return bool(len(positions) and min(np.hypot(positions[:,0]-y,positions[:,1]-x))<1.5)


def exposure(row,ra,dec,source):
    src_ra=ra+float(source['east_arcsec'])/(3600*np.cos(np.deg2rad(dec)))
    src_dec=dec+float(source['north_arcsec'])/3600
    _,im,flag,w,left,down=cutout(row['url'],ra,dec,[])
    def xy(a,d):
        x,y=w.all_world2pix(a,d,0)
        return float(x)-left,float(y)-down
    sx,sy=xy(src_ra,src_dec)
    valid=((flag & ~((1<<18)|(1<<24)))==0)&np.isfinite(im)
    smooth=gaussian_filter(im,.65)
    resid=smooth-median_filter(smooth,11)
    med=np.median(resid[valid])
    threshold=med+6*1.4826*np.median(abs(resid[valid]-med))
    old=peaks(resid,valid,threshold)
    if not nearby(old,sx,sy):return {},'assumed source not recovered'
    iy,ix=np.rint([sy,sx]).astype(int)
    yy,xx=np.indices(im.shape);rad=np.hypot(yy-sy,xx-sx)
    if not np.all(valid[rad<=1.8]):return {},'masked source core'
    bg=np.median(im[valid&(rad>=3.5)&(rad<=5.5)])
    stamp=np.where(valid&(rad<=3.3),np.maximum(im-bg,0.),0.)[iy-4:iy+5,ix-4:ix+5]
    yg,xg=np.indices(stamp.shape)
    core=stamp*(np.hypot(yg-4,xg-4)<=2)
    if core.sum()<=0:return {},'empty source core'
    cx=np.average(xg,weights=core);cy=np.average(yg,weights=core)
    result={}
    for radius in RADII:
        for angle in range(NANGLE):
            phi=angle*2*np.pi/NANGLE
            x,y=xy(src_ra+radius*np.sin(phi)/(3600*np.cos(np.deg2rad(src_dec))),
                   src_dec+radius*np.cos(phi)/3600)
            j,i=np.rint([x,y]).astype(int)
            blocked=nearby(old,x,y) or not np.all(valid[np.hypot(yy-y,xx-x)<=1.8])
            addition=np.zeros_like(im)
            shifted=shift(stamp,(y-i+4-cy,x-j+4-cx),order=1,mode='constant',prefilter=False)
            if shifted.sum()>0:shifted*=stamp.sum()/shifted.sum()
            addition[i-4:i+5,j-4:j+5]=shifted
            filtered=gaussian_filter(addition,.65)
            for delta in DELTAS:
                r=resid+filtered*10**(-.4*delta)
                result[(radius,delta,angle)]=not blocked and nearby(peaks(r,valid,threshold),x,y)
    return result,'ok'


if __name__=='__main__':
    CACHE.mkdir(exist_ok=True,parents=True)
    with (HERE/'other/euclid_scene_checks.csv').open() as f:frames=list(csv.DictReader(f))
    with (HERE/'other/euclid_scene_sources.csv').open() as f:sources=list(csv.DictReader(f))
    extra=json.loads((HERE/'large-source-images.json').read_text())
    frames+=extra['euclid_queries'];sources+=extra['euclid_sources']
    with (HERE/'kmt/catalog-2024.tsv').open() as f:
        positions={r['Event']:coords(r['RA'],r['Dec']) for r in csv.DictReader(f,delimiter='\t')}
    output=[];receipts=[]
    for event in ['KMT-2024-BLG-0456','KMT-2024-BLG-3278']:
        opt=min((s for s in sources if s['event']==event and int(s['independent_dithers'])>=3),
                key=lambda s:float(s['radius_arcsec']))
        detections=[]
        for row in frames:
            if row['event']!=event or row['status']!='ok':continue
            recovered,status=exposure(row,*positions[event],opt)
            receipts.append(dict(event=event,url=row['url'],mjd=row['mjd'],status=status))
            if status=='ok':detections.append(recovered)
        if len(detections)<2:raise ValueError('too few frames containing the assumed source')
        minimum=max(2,math.ceil(len(detections)/2))
        for rad in RADII:
            for delta in DELTAS:
                hits=sum(sum(int(r[rad,delta,i]) for r in detections)>=minimum for i in range(NANGLE))
                output.append(dict(event=event,source_offset_arcsec=opt['radius_arcsec'],separation_arcsec=rad,
                                   delta_vis_mag=delta,qualifying_frames=len(detections),minimum_frame_detections=minimum,
                                   position_angles=NANGLE,recovered_angles=hits,fraction=round(hits/NANGLE,4)))
    with (HERE/'tables/conditional-image-contrast.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(output[0]),lineterminator='\n');w.writeheader();w.writerows(output)
    (HERE/'tables/image-contrast-frames.json').write_text(json.dumps(receipts,indent=2)+'\n')
