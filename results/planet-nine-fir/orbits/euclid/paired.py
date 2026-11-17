"""Run native VIS injections on the identical DES off-grid trials and join the outcome."""
import argparse
import csv
import gzip
import hashlib
import json
import numpy as np
from astropy.time import Time
from scipy.ndimage import maximum_filter,map_coordinates
from scipy.constants import astronomical_unit

from dense import BASE,CACHE,HERE,OUT,R,design,kernel,p
from mosaic import xy
from vis_pixels import TREF,SC
from offgrid import catalog,scores,wilson

MODELS=[('mass5_R2.71_p0.10',2.71,.1),('mass5_R2.94_p0.30',2.94,.3),
        ('mass10_R3.70_p0.30',3.7,.3)]

def main():
    ap=argparse.ArgumentParser()
    from pathlib import Path
    ap.add_argument("des_draws",type=Path)
    args=ap.parse_args()
    rows=list(csv.DictReader(gzip.open(args.des_draws,'rt')))
    vec=lambda keys:np.array([[float(row[k]) for k in keys] for row in rows])
    r=vec(['x_au','y_au','z_au']);v=vec(['vx_au_d','vy_au_d','vz_au_d'])
    t0=vec(['reference_emission_tt_days_since_j2000'])[:,0]
    states=(r,v,t0,{})
    t=TREF.tt.jd-2451545.; obs=SC(TREF.tdb.jd-2451545.)
    ra,de=p.sky(states,t,obs)
    dt=t-t0-p.norm(r-obs)/p.C
    for i in range(3):
        pos=p.kepler(r,v,dt)
        delta=p.norm(pos-obs)
        dt=t-t0-delta/p.C
    heli_earth,bary_earth=p.observers(t0+dt)
    heli=p.norm(pos-(bary_earth-heli_earth))
    solar=json.loads((HERE/'solar_vis.json').read_text())['sun_VIS_AB_at_1au']
    mags=np.array([solar-2.5*np.log10(albedo*(radius*6371e3/astronomical_unit/(heli*delta))**2)
                   for key,radius,albedo in MODELS])
    names,ts,bank,_=design()
    ks=np.array([[np.outer(kernel(dy),kernel(dx)) for dx,dy in motion] for motion in bank])
    g=json.loads((BASE/'geometry.json').read_text())
    D=[np.load(BASE/(n+'_data.npy'),mmap_mode='r') for n in names]
    W=[np.load(BASE/(n+'_weight.npy'),mmap_mode='r') for n in names]
    field=np.load(OUT/'eligible.npy',mmap_mode='r')
    hull=np.load(OUT/'support.npy',mmap_mode='r')
    z=np.load(OUT/'score.npy',mmap_mode='r')
    threshold=np.maximum(8.,maximum_filter(z,size=25)+.01)
    psf=np.load(CACHE/'measured_psf.npy')
    receipt=list(csv.DictReader((HERE/'mosaic_wcs.csv').open()))
    counts=np.array([3.48*np.median([float(q['ADU_at_VIS25']) for q in receipt if q['quadrant'].startswith(n)]) for n in names])
    xs,ys=xy(ra,de);xs-=g['x0'];ys-=g['y0']
    ix,iy=np.rint(xs).astype(int),np.rint(ys).astype(int)
    safe=(ix>=R)&(ix<g['width']-R)&(iy>=R)&(iy<g['height']-R)
    i=np.clip(ix,0,g['width']-1);j=np.clip(iy,0,g['height']-1)
    sep,_=catalog().query(np.column_stack((xs,ys)))
    eligible=safe&(hull[j,i]!=0)&(field[j,i]!=0)&(sep>15)
    drift=[]
    for mjd in ts:
        tt=Time(mjd,format='mjd',scale='utc')
        a,b=p.sky(states,tt.tt.jd-2451545,SC(tt.tdb.jd-2451545))
        x,y=xy(a,b)
        drift.append((x-g['x0']-ix,y-g['y0']-iy))
    dx,dy=np.meshgrid(np.arange(-R,R+1),np.arange(-R,R+1))
    ids=np.flatnonzero(eligible)
    vis=np.zeros((len(r),3),dtype=bool)
    for offset in range(0,len(ids),48):
        sel=ids[offset:offset+48]
        loc=(iy[sel,None,None]+dy,ix[sel,None,None]+dx)
        data=np.array([q[loc] for q in D]).transpose(1,0,2,3)
        weight=np.array([q[loc] for q in W]).transpose(1,0,2,3)
        src=[]
        for f,(x,y) in enumerate(drift):
            coords=[dy[None,:,:]-y[sel,None,None]+psf.shape[1]//2,
                    dx[None,:,:]-x[sel,None,None]+psf.shape[2]//2]
            src.append(map_coordinates(psf[f],coords,order=1,cval=0,prefilter=False))
        src=np.array(src).transpose(1,0,2,3)
        vis[sel],_=scores(data,weight,src,ks,counts,20260316+offset,mags[:,sel],threshold[iy[sel],ix[sel]])
    des=np.array([[bool(int(row['hit_'+model[0]])) for model in MODELS] for row in rows])
    union=des|vis
    with (HERE/'paired_trials.csv').open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['experiment','draw','euclid_ra_deg','euclid_dec_deg','euclid_heliocentric_au','euclid_range_au','euclid_eligible']+[label+'_'+m[0] for m in MODELS for label in ('VIS_AB','DES_found','Euclid_found','either_found')])
        for k,row in enumerate(rows):
            values=[]
            for im in range(3): values.extend([round(float(mags[im,k]),5),int(des[k,im]),int(vis[k,im]),int(union[k,im])])
            w.writerow([row['experiment'],row['draw'],ra[k],de[k],heli[k],delta[k],int(eligible[k]),*values])
    summary=dict(des_input_sha256=hashlib.sha256(args.des_draws.read_bytes()).hexdigest(),cases={})
    exp_names=np.array([q['experiment'] for q in rows])
    for group in np.unique(exp_names):
        keep=exp_names==group
        summary['cases'][group]={}
        for im,(key,*_) in enumerate(MODELS):
            both=des[:,im]&vis[:,im]
            k=int(union[keep,im].sum());n=int(keep.sum())
            summary['cases'][group][key]=dict(n=n,des=int(des[keep,im].sum()),euclid=int(vis[keep,im].sum()),
                both=int(both[keep].sum()),either=k,neither=n-k,wilson95=wilson(k,n))
    (HERE/'paired_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary)


if __name__=='__main__': main()
