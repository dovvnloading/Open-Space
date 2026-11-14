"""Fit models directly to target photometry and evaluate held-out time blocks.

The independent-clock and four-dip blend models add eclipse deficits in linear
flux. Each model includes an independent out-of-eclipse level per OGLE survey.
"""
import argparse
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
PERIOD = 1.1696719
EPOCH = 2460000.0369
TIME_ORIGIN = 2456000.0
INPUTS = {
    'ogle3.dat': ['https://www.astrouw.edu.pl/ogle/ogle3/OIII-CVS/smc/ecl/phot/I/OGLE-SMC-ECL-6090.dat',
                  'd57614117ef140adf1413baab42f92182b9c39752eed39706b6fc040e203a897'],
    'ogle4.dat': ['https://www.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/phot/I/OGLE-SMC-ECL-06090.dat',
                  '429caaf74f2086f3a8f53c71a91bfb3c07049dbc2e81122c06646387c8fa21a0'],
}
def get_data():
    pieces=[]
    for survey,(name,(_,sha)) in enumerate(INPUTS.items()):
        path=HERE/name
        assert hashlib.sha256(path.read_bytes()).hexdigest()==sha
        t,m,e=np.loadtxt(path).T
        keep=np.isfinite(t+m+e)&(e>0)&(e<=.10)
        t,m,e=t[keep],m[keep],e[keep]
        if survey==0:t=t+2450000
        flux=10**(-.4*(m-18.23))
        err=np.hypot(e,.012)*(.4*np.log(10))*flux
        pieces.append(np.column_stack((t,(t-TIME_ORIGIN)/365.25,((t-EPOCH)/PERIOD)%1,
                                      np.full(len(t),survey),flux,err)))
    return np.concatenate(pieces)
def wrapped(x):
    return (x+.5)%1-.5

def solve_model(p,D,kind,jac=False):
    _,year,phase,survey,_,_=D.T
    y=np.where(survey==0,p[0],p[1])
    if jac:
        J=np.zeros((len(D),len(p)))
        J[:,0]=survey==0
        J[:,1]=survey==1
    def dip(center,width,depth,ic,iw,idp,slope,extra=None):
        u=wrapped(phase-center)/width
        kernel=np.exp(-.5*u*u)
        if jac:
            dc=-depth*kernel*u/width
            J[:,ic]+=dc
            J[:,slope]+=year*dc
            J[:,iw]-=depth*kernel*u*u/width
            J[:,idp]-=kernel
            if extra is not None:
                for key,mult in extra:
                    J[:,key]+=dc*mult
        y[:]-=depth*kernel
    if kind=='apsidal':
        ang=2*np.pi*year/p[12]
        co,si=np.cos(ang),np.sin(ang)
        move=p[10]*co+p[11]*si
        dper=(p[10]*si-p[11]*co)*2*np.pi*year/p[12]**2
        for sign,i,j,k in [(-1,2,4,6),(1,3,5,7)]:
            dip(p[i]+p[8]*year+sign*move,p[j],p[k],i,j,k,8,
                [(10,sign*co),(11,sign*si),(12,sign*dper)])
    elif kind=='blend':
        for i,j,k,rate in [(2,4,6,8),(3,5,7,8),(10,12,14,16),(11,13,15,16)]:
            dip(p[i]+p[rate]*year,p[j],p[k],i,j,k,rate)
    else:
        for i,j,k,rate in [(2,4,6,8),(3,5,7,8 if kind=='static' else 10)]:
            dip(p[i]+p[rate]*year,p[j],p[k],i,j,k,rate)
    # A common static quadrature variation accommodates ellipsoidal background in all fits.
    h=np.cos(4*np.pi*phase)
    y+=p[9]*h
    if jac:J[:,9]=h
    return J if jac else y

def setup(kind):
    base=[1,1,.01,.51,.024,.023,.32,.31,0,.0]
    low=[.90,.90,-.2,.28,.006,.006,.01,.01,-.055,-.06]
    high=[1.10,1.10,.23,.75,.095,.095,.75,.75,.055,.06]
    if kind=='apsidal':
        base += [-.020,0,32]
        low += [-.13,-.13,6]
        high += [.13,.13,150]
    elif kind=='two_clocks':
        base+=[.002]; low+=[-.055]; high+=[.055]
    elif kind=='blend':
        # Eclipse locations, widths and attenuated depths fully independent between binaries.
        # Requiring one primary near phase 0 and one secondary near phase .5 per binary tests close-period blends.
        base += [.01,.52,.021,.023,.12,.12,.008]
        low += [-.36,.12,.006,.006,0,0,-.11]
        high += [.38,.89,.095,.095,.75,.75,.11]
        low[6]=low[7]=.0
        base[6]=base[7]=.22
    return np.array(base), np.array(low), np.array(high)

def starts(kind,full):
    base,lo,hi=setup(kind)
    if kind=='apsidal':
        for per in [12,18,25,32,42,60,100]:
            p=base.copy();p[12]=per;yield p
    elif kind=='blend':
        rng=np.random.default_rng(636412)
        yield base.copy()
        for i in range(95 if full else 35):
            p=base.copy()
            if i%4==0:
                # includes one primary eclipsing clock and a separate secondary clock; counterparts may vanish.
                p[7]=.015;p[6]=.34;p[14]=.02;p[15]=.30
                p[3]=rng.uniform(.35,.71);p[10]=rng.uniform(-.19,.17)
            elif i%4==1:
                p[2:4]+=rng.normal(0,.035,2);p[10:12]+=rng.normal(0,.055,2)
                p[[6,7,14,15]]=rng.uniform(.07,.27,4)
            elif i%4==2:
                p[10:12]+=rng.normal(0,.10,2)
                p[[14,15]]=rng.uniform(.02,.11,2)
            else:
                f=rng.uniform(.18,.8,2)
                p[6:8]=.34*f;p[14:16]=.34*(1-f)
                p[2:4]+=rng.normal(0,.025,2);p[10:12]+=rng.normal(0,.025,2)
            p[8]=rng.uniform(-.010,.010);p[16]=rng.uniform(-.018,.018) if i%3 else rng.uniform(-.085,.085)
            yield np.clip(p,lo+1e-6,hi-1e-6)
    else:
        yield base
        # Independent near-catalog clocks can drift by up to ~1.2 cycles over the observations.
        for s in [-.018,-.008,.004,.012,.022]:
            p=base.copy();p[8]=s
            if kind=='two_clocks':p[10]=-s/2
            yield p

def fit_job(task):
    fold,kind,full=task
    D=get_data();codes=fold_id(D[:,0])
    train=np.ones(len(D),bool) if fold=='all' else codes!=fold
    T=D[train];base,lo,hi=setup(kind)
    def fun(p):return (solve_model(p,T,kind)-T[:,4])/T[:,5]
    def deriv(p):return solve_model(p,T,kind,True)/T[:,5,None]
    best=None; costs=[]
    for p0 in starts(kind,full):
        f=least_squares(fun,p0,jac=deriv,bounds=(lo,hi),loss='soft_l1',f_scale=1.8,
                        ftol=3e-7,xtol=3e-7,gtol=3e-7,max_nfev=400, x_scale='jac')
        costs.append(float(f.cost))
        if best is None or f.cost<best.cost:best=f
    r=(solve_model(best.x,D,kind)-D[:,4])/D[:,5]
    metrics={}
    for name,mask in [('train',train),('test',~train)]:
        if not mask.any():continue
        err=(solve_model(best.x,D[mask],kind)-D[mask,4])
        metrics[name]={'n':int(mask.sum()),'chi2':float(np.sum(r[mask]**2)),
                       'soft_deviance':float(2*1.8**2*np.sum(np.sqrt(1+(r[mask]/1.8)**2)-1)),
                       'rms_flux':float(np.sqrt(np.mean(err**2)))}
    params=best.x.tolist()
    print(fold,kind,round(2*best.cost,1),metrics.get('test',{}).get('soft_deviance'),flush=True)
    return dict(fold=fold,model=kind,params=params,metrics=metrics,
                starts=len(costs),five_best_training_costs=sorted(costs)[:5],success=bool(best.success))
def fold_id(t):
    bounds=[2453500,2455200,2457000,2458300,2459700]
    return np.searchsorted(bounds,t,side='right')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--model',choices=['all','blend','basic'],default='all')
    parser.add_argument('--output',type=Path,default=HERE/'fits.json')
    arg=parser.parse_args()
    models={'all':['static','two_clocks','apsidal','blend'],'basic':['static','two_clocks','apsidal'],'blend':['blend']}[arg.model]
    tasks=[(fold,kind,fold=='all') for fold in ['all',5,0,2,1,4,3] for kind in models]
    with ProcessPoolExecutor(max_workers=arg.workers) as pool: res=list(pool.map(fit_job,tasks))
    arg.output.write_text(json.dumps({'input_urls_and_sha256':INPUTS,'fitted':res},indent=2)+'\n')
if __name__=='__main__':main()
