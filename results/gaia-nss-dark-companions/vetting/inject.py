"""Empirical recovery of late-dwarf secondary absorption in two LAMOST blue coadds."""
import argparse
import csv
import hashlib
import json
import urllib.request
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d, maximum_filter1d, percentile_filter
from spectra import C, HERE, IDS, read_header, spectrum

BASE='https://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/'
FILES={'wave':'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits',
       'kstar':'PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte05400-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'}
FRACTIONS=[0,.01,.02,.026,.036,.05,.1]
# Offsets much wider than one resolution element sample unrelated structures in the real residuals.
OFFSETS=np.r_[np.arange(-1800,-199,60),0,np.arange(240,1801,60)]


def load_template(cache,vsini,arm):
    arrays=[]; meta={}
    for name,suffix in FILES.items():
        dest=cache/(name+'.fits')
        if not dest.exists():dest.write_bytes(urllib.request.urlopen(BASE+suffix,timeout=120).read())
        with dest.open('rb') as f:
            h=read_header(f)
            arrays.append(np.fromfile(f,dtype=('>f8' if h['BITPIX']=='-64' else '>f4'),count=int(h['NAXIS1'])))
        meta[name]=dict(url=BASE+suffix,sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
    wave,flux=arrays
    lo,hi=(4850,5450) if arm=='B' else (6250,6950)
    x=np.exp(np.arange(np.log(lo),np.log(hi),2/C))
    y=np.interp(x,wave,flux)
    # Linear limb darkening epsilon=0.6; x samples 2 km/s.
    v=np.arange(-int(vsini/2),int(vsini/2)+1)*2/vsini
    k=2*.4*np.sqrt(1-v*v)+.5*np.pi*.6*(1-v*v)
    y=np.convolve(y,k/k.sum(),mode='same')
    y=gaussian_filter1d(y,C/7500/2.355/2)
    continuum=gaussian_filter1d(percentile_filter(y,85,size=351),65)
    return x,1-y/continuum,meta


def residualize(v,basis):
    return v-basis@(basis.T@v)


def run(cache,pair=IDS,vel_grid=OFFSETS,arm="B",weak_only=False):
    from scipy.optimize import minimize_scalar
    from spectra import URL
    products=[]
    for obsid in pair:
        path=cache/(obsid+'.fits.gz')
        if not path.exists():
            path.write_bytes(urllib.request.urlopen(URL+obsid,timeout=90).read())
        products.append(dict(obsid=obsid,url=URL+obsid,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    (x1,y1,g1),(x2,y2,g2)=[spectrum(cache/(oid+'.fits.gz'),arm) for oid in pair]
    lo,hi=(4950,5340) if arm=='B' else (6400,6800)
    use=(x2>=lo)&(x2<=hi)
    if arm=='R': use &= ((x2<6510)|(x2>6615))
    core=use&g2
    for u in np.arange(0,61,2):core &= np.interp(x2/(1+u/C),x1,g1.astype(float))>.999
    def mismatch(u):
        t=np.interp(x2[core]/(1+u/C),x1,y1)
        return -np.corrcoef(t,y2[core])[0,1]
    dv=minimize_scalar(mismatch,bounds=(0,60),method='bounded').x
    mask=use&g2&(np.interp(x2/(1+dv/C),x1,g1.astype(float))>.999)
    if weak_only:
        # Reject primary lines deeper than 5% within approximately one resolution element
        # of a pixel in either spectrum. Leave the same mask fixed at all trial velocities.
        m1=maximum_filter1d(np.abs(y1),size=11)
        m2=maximum_filter1d(np.abs(y2),size=11)
        mask &= (m2<.05)&(np.interp(x2/(1+dv/C),x1,m1)<.05)
    x=x2[mask]; low_x=x/(1+dv/C)
    low=np.interp(low_x,x1,y1);high=y2[mask]
    # Do not force a constant continuum, line-spread function or wavelength offset across the blue arm.
    edges=np.linspace(lo,hi,11)
    blocks=[]
    for a,b in zip(edges,edges[1:]):
        ix=(x>=a)&(x<b);w=x[ix]
        if len(w)<6:continue
        lsignal=np.interp(low_x[ix],x1,y1)
        d1=np.gradient(lsignal,w)
        d2=np.gradient(d1,w)
        design=np.column_stack((np.ones(len(w)),(w-w.mean())/w.std(),lsignal,d1,d2))
        q,_=np.linalg.qr(design)
        blocks.append((ix,q))
    def remove_nuisance(y):
        out=np.zeros_like(y)
        for ix,q in blocks:out[ix]=residualize(y[ix],q)
        return out
    base=remove_nuisance(high-low)
    source='3378588057203660160'
    rows=[r for r in csv.DictReader((HERE.parent/'lamost/exposures.csv').open()) if r['source_id']==source]
    prim=[np.median([float(r['rv_kms']) for r in rows if r['obsid']==o]) for o in pair]
    # Systemic velocity and representative masses from the joint orbit.
    gamma=32.564
    sec=[gamma-(1.95/1.77)*(v-gamma) for v in prim]
    scenarios=[]; nullrows=[]
    for rot in [6.,40.]:
        grid,model,meta=load_template(cache,rot,arm)
        for v0,v1 in [(0,0),(40,40),(0,80),(80,120)]:
            def twins(wave,velocity,split):
                return (np.interp(wave/(1+(velocity-split)/C),grid,model)+
                        np.interp(wave/(1+(velocity+split)/C),grid,model))/2
            sh=[]
            for off in vel_grid:
                cool_low=twins(low_x,sec[0]+off,v0)
                cool_high=twins(x,sec[1]+off,v1)
                sh.append(remove_nuisance(cool_high-cool_low))
            sh=np.array(sh)
            fn=np.sum(base*sh,axis=1)/np.sum(sh*sh,axis=1)
            i0=np.flatnonzero(vel_grid==0).item()
            null=fn[np.abs(vel_grid)>=200]
            q95=float(np.quantile(null,.95))
            # In flux units, mixing a secondary at fraction f in both actual spectra yields
            # (1-f)*base + f*shape after differencing and nuisance projection.
            powers={}
            for frac in FRACTIONS:
                injected=(1-frac)*base+frac*sh
                measured=np.sum(injected*sh,axis=1)/np.sum(sh*sh,axis=1)
                controls=measured[np.abs(vel_grid)>=200]
                powers[str(frac)]=float(np.mean([v>np.quantile(np.delete(null,i),.95)
                                                  for i,v in enumerate(controls)]))
            local=np.abs(vel_grid)<=150
            best_idx=np.argmax(np.where(local,fn,-np.inf))
            scenarios.append(dict(local_max_offset_kms=int(vel_grid[best_idx]),local_max_fraction=float(fn[best_idx]),
                 vsini=rot,inner_half_separations_kms=[v0,v1],n_controls=len(null),
                 target_fitted_fraction=float(fn[i0]),control_p95=q95,
                 target_exceeds_control_p95=bool(fn[i0]>q95),injected_fraction_recovery=powers))
            for off,f in zip(vel_grid,fn):
                nullrows.append(dict(vsini=rot,v_low=v0,v_high=v1,offset_kms=int(off),fitted_fraction=f))
    if not weak_only:
        with (HERE/(('red_' if arm=='R' else '')+('injection_control_scores.csv' if pair==IDS else 'control_scores_'+pair[0]+'.csv'))).open('w',newline='') as file:
            w=csv.DictWriter(file,list(nullrows[0]),lineterminator='\n');w.writeheader();w.writerows(nullrows)
    return dict(source_id=source,coadds=products,measured_shift_kms=dv,model_files=meta,resolving_power=7500,arm=arm,wavelength_A=[lo,hi],
            matched_pixels=len(x),weak_primary_mask=weak_only,rv_primary_kms=prim,rv_secondary_barycentre_kms=sec,
            scenarios=scenarios)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',required=True,type=Path)
    a=p.parse_args()
    (HERE/'injections.json').write_text(json.dumps(run(a.cache,IDS,np.unique(np.r_[OFFSETS,np.arange(-150,151,3)])),indent=2)+'\n')
    checks={}
    for pair in [('653004120','788407188'),('638204120','786707188'),('895707188','799107188')]:
        checks[pair[0]]=run(a.cache,pair,np.unique(np.r_[OFFSETS,np.arange(-150,151,3)]))
    (HERE/'replicated_injections.json').write_text(json.dumps(checks,indent=2)+'\n')
    red={}
    for pair in [IDS,('653004120','788407188'),('638204120','786707188'),('895707188','799107188')]:
        red[pair[0]]=run(a.cache,pair,np.unique(np.r_[OFFSETS,np.arange(-150,151,3)]),arm='R')
    (HERE/'red_injections.json').write_text(json.dumps(red,indent=2)+'\n')
    weak={}
    for arm in ['B','R']:
        for pair in [IDS,('653004120','788407188'),('638204120','786707188'),('895707188','799107188')]:
            weak[arm+pair[0]]=run(a.cache,pair,OFFSETS,arm=arm,weak_only=True)
    (HERE/'weak_line_injections.json').write_text(json.dumps(weak,indent=2)+'\n')
