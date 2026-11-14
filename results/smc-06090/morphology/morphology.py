"""Measure both flanks of both eclipses in each OGLE time block, profiling skewness.

Run this file next to the two original source tables in inputs/ or supply --inputs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

P=1.1696719
T0=2460000.0369
HERE=Path(__file__).resolve().parent
FILES={
 'ogle3.dat':('d57614117ef140adf1413baab42f92182b9c39752eed39706b6fc040e203a897','https://www.astrouw.edu.pl/ogle/ogle3/OIII-CVS/smc/ecl/phot/I/OGLE-SMC-ECL-6090.dat'),
 'ogle4.dat':('429caaf74f2086f3a8f53c71a91bfb3c07049dbc2e81122c06646387c8fa21a0','https://www.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/phot/I/OGLE-SMC-ECL-06090.dat')}
BINS=[('ogle3.dat',0,2453500),('ogle3.dat',2453500,np.inf),
      ('ogle4.dat',0,2457000),('ogle4.dat',2457000,2458300),
      ('ogle4.dat',2458300,2459700),('ogle4.dat',2459700,np.inf)]
FIELDS=['center_separation','halfdepth_midpoint_separation','ingress_separation','egress_separation',
        'primary_fwhm_cycles','secondary_fwhm_cycles','primary_depth_mag','secondary_depth_mag',
        'primary_halfwidth_imbalance','secondary_halfwidth_imbalance']

def wrap(x):
    return (x+.5)%1-.5

def model(z,phase):
    out=np.full_like(phase,z[0])
    for mag,center,left,right in (z[1:5],z[5:9]):
        u=wrap(phase-center)
        out+=mag*np.exp(-.5*(u/np.where(u<0,left,right))**2)
    return out

def metrics(z):
    _,depth1,c1,l1,r1,depth2,c2,l2,r2=z
    f=np.sqrt(2*np.log(2))
    left=(c2-f*l2)-(c1-f*l1)
    right=(c2+f*r2)-(c1+f*r1)
    return np.array([c2-c1,(left+right)/2,left,right,f*(l1+r1),f*(l2+r2),
                     depth1,depth2,(r1-l1)/(r1+l1),(r2-l2)/(r2+l2)])

def fit(p,m,e,initial=None):
    b=np.quantile(m,.25)
    if initial is None: initial=[b,.4,.015,.025,.025,.4,.51,.025,.025]
    lower=[b-.3,.005,-.14,.004,.004,.005,.35,.004,.004]
    upper=[b+.3,1.5,.14,.10,.10,1.5,.65,.10,.10]
    result=least_squares(lambda z:(model(z,p)-m)/np.hypot(e,.012),initial,
                        bounds=(lower,upper),loss='soft_l1',f_scale=1.8,max_nfev=180)
    return result.x

def run(inputs,draws):
    rng=np.random.default_rng(92763);out=[];sim=[];point=[]
    for name,lo,hi in BINS:
        path=inputs/name
        assert hashlib.sha256(path.read_bytes()).hexdigest()==FILES[name][0]
        t,m,e=np.loadtxt(path,usecols=(0,1,2)).T
        if name=='ogle3.dat':t+=2450000
        ok=np.isfinite(t+m+e)&(e>0)&(e<=.10)&(t>=lo)&(t<hi)
        t,m,e=t[ok],m[ok],e[ok];ph=wrap((t-T0)/P)
        z=fit(ph,m,e);q=metrics(z)
        _,ix=np.unique(np.floor(t),return_inverse=True)
        groups=[np.flatnonzero(ix==k) for k in range(1+ix.max())]
        samples=[]
        bounds=0
        for _ in range(draws):
            ind=np.concatenate([groups[k] for k in rng.integers(0,len(groups),len(groups))])
            params=fit(ph[ind],m[ind],e[ind],z)
            w=params[[3,4,7,8]]
            bounds+=bool((w<.00401).any() or (w>.09999).any())
            samples.append(metrics(params))
        a=np.array(samples)
        point.append(q);sim.append(a)
        out.append(dict(observations=len(t),hjd_range=[t.min(),t.max()],params=z.tolist(),
                        values=q.tolist(),ci95=np.percentile(a,[2.5,97.5],axis=0).T.tolist(),
                        draws_with_width_at_bound=bounds))
    diffs={}
    for key,start,end in [('OGLE-III',0,1),('OGLE-IV',2,5)]:
        diffs[key]=dict(delta=(point[end]-point[start]).tolist(),
                         ci95=np.percentile(sim[end]-sim[start],[2.5,97.5],axis=0).T.tolist())
    return dict(models='Separate ingress and egress Gaussian half-widths for each eclipse; seven-parameter symmetric fit extended by two widths.',inputs=FILES,period=P,hjd_epoch=T0,draws=draws,seed=92763,fields=FIELDS,epochs=out,changes=diffs)

if __name__=='__main__':
    cli=argparse.ArgumentParser();cli.add_argument('--inputs',type=Path,default=HERE.parent/'joint'/'inputs')
    cli.add_argument('--draws',type=int,default=600)
    cli.add_argument('--out',type=Path,default=HERE/'morphology.json')
    args=cli.parse_args()
    r=run(args.inputs,args.draws);args.out.write_text(json.dumps(r,indent=2)+'\n')
    for k,v in r['changes'].items():print(k,np.round(v['delta'][:4],5),np.round(v['ci95'][:4],5))
