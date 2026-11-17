"""Continuously optimize single-lens and asymmetric exponential models on short leads.

Fits use all valid epochs, independent stream zero points, and a robust
off-peak error inflation factor. Screen on disk with run_grid.py first.
"""
import argparse
import csv
import json
import multiprocessing as mp
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import fit
from run_grid import inputs, read_curve


def local_fit(job, seed):
    vb=fit.engine()
    name,_,anchor,guess=job
    t,f,e,g,diff=read_curve(job)
    t,f,e,g=fit.clean(t,f,e,g)
    if diff:
        e=fit.reported_error_floor(e,g)
    origin=seed['fs_t0']
    t=t-origin
    ng=max(g)+1
    # Estimate a multiplicative scale from the central 68% of standardized
    # residuals off peak (never deflate). Retain every epoch in the fit.
    cap=np.full(ng,1e28)
    for j in range(ng):
        inside=g==j
        quiet=inside&(abs(t)>max(5,min(2*abs(guess),100)) if np.isfinite(guess) else abs(t)>5)
        if quiet.sum()>=10:
            scale=np.median(abs((f[quiet]-np.median(f[quiet]))/e[quiet]))/0.67448975
            e[inside]*=max(1.,scale)
        if not diff:
            cap[j]=max(1e-12,1.15*np.median(f[inside]))
    w=1/e**2
    sums=np.bincount(g,weights=w)
    mean=np.bincount(g,weights=w*f)/sums
    target=f-mean[g]
    null=np.sum(w*target**2)

    def pred(x,kind):
        dt,lte,lu=x[:3]
        te,u=np.exp([lte,lu])
        tau=(t-dt)/te
        if kind=='flare':
            return np.exp(np.clip(np.where(t<dt,(t-dt)/te,-(t-dt)/u),-745,0))
        sep=np.hypot(tau,u)
        magn=fit.lens.magnification(sep**2)
        if kind=='fs':
            rho=np.exp(x[3])
            close=sep<max(6*rho,2)
            if np.any(close):
                # At greater separations the point-source limit is used.
                magn[close]=vb.ESPLLightCurve([lu,lte,dt,x[3]],t[close].tolist())[0]
        return magn-1

    def residual(x,kind,details=False):
        a=pred(x,kind)
        avg=np.bincount(g,weights=w*a)/sums
        ax=a-avg[g]
        xx=np.bincount(g,weights=w*ax**2)
        xy=np.bincount(g,weights=w*ax*target)
        fs=np.clip(xy/np.maximum(1e-300,xx),0.,np.inf if kind=='flare' else cap)
        rr=(target-fs[g]*ax)/e
        if details: return rr,fs,mean-fs*avg
        return np.sum(rr**2)

    # Optimize relative to a reference HJD for stable differences.
    bounds=[(-max(4,min(30,abs(guess) if np.isfinite(guess) else 4)),
              max(4,min(30,abs(guess) if np.isfinite(guess) else 4))),
            (np.log(.007),np.log(1000)),(np.log(.002),np.log(12))]
    output={'event':name,'n':len(t),'streams':int(ng),'null_chi2':round(null,3)}
    use=seed
    best={}
    for kind in ['ps','fs','flare']:
        if kind=='flare':
            starts=[[0,np.log(.08),np.log(.2)],[0,np.log(.3),np.log(.06)]]
            limits=bounds[:2]+[(np.log(.003),np.log(200))]
        else:
            starts=[]
            for fam in ('ps','fs'):
                starts.append([use[fam+'_t0']-origin,np.log(use[fam+'_te']),np.log(use[fam+'_u0'])])
            if np.isfinite(guess) and guess>0 and np.isfinite(anchor):
                starts.append([anchor-origin,np.log(guess),np.log(.2)])
            if kind=='fs':
                starts=[x+[np.log(r)] for x,r in zip(starts,[max(.1,use['fs_rho']),use['fs_rho'],1.])]
                starts.append(list(best['ps'].x)+[np.log(.08)])
                starts.append([0.,np.log(.22),np.log(.45),np.log(2.0)])
            limits=bounds+([(np.log(.01),np.log(20))] if kind=='fs' else [])
        tries=[]
        for start in starts:
            start=np.clip(start,[i[0] for i in limits],[i[1] for i in limits])
            # The simplex must be able to move a central epoch initially at zero.
            simplex=np.tile(start,(len(start)+1,1))
            for i in range(len(start)):
                step=.04 if i==0 else .18
                simplex[i+1,i]+=step if start[i]+step<limits[i][1] else -step
            sol=minimize(residual,start,args=(kind,),method='Nelder-Mead',bounds=limits,
                         options={'maxiter':300,'xatol':0.0004,'fatol':0.008,'initial_simplex':simplex})
            tries.append(sol)
        b=best[kind]=min(tries,key=lambda s:s.fun)
        output[kind+'_chi2']=round(b.fun,3)
        output[kind+'_t0']=round(b.x[0]+origin,6)
        output[kind+('_rise_d' if kind=='flare' else '_te')]=round(np.exp(b.x[1]),6)
        output[kind+('_fall_d' if kind=='flare' else '_u0')]=round(np.exp(b.x[2]),6)
        if kind=='fs':
            output['rho']=round(np.exp(b.x[3]),6)
            out,fs,base=residual(b.x,'fs',True)
            window=abs(t-b.x[0])<max(.18,1.5*np.exp(b.x[1])*max(1.,np.exp(b.x[3])))
            output['peak_n']=int(window.sum());output['peak_streams']=int(np.unique(g[window]).size)
            output['peak_chi2']=round(np.sum(out[window]**2),3)
            output['base_chi2_per_epoch']=round(np.mean(out[~window]**2),4) if np.any(~window) else None
    return output


def worker(arg):
    try:
        return local_fit(*arg)
    except (ValueError,IndexError,ArithmeticError,OSError) as ex:
        return {'event':arg[0][0],'error':str(ex)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('names',nargs='+');p.add_argument('--cache',type=Path,default=Path.cwd())
    p.add_argument('--grid',type=Path,default=Path('staging/grid'))
    p.add_argument('--out',type=Path,default=Path('staging/refined'))
    p.add_argument('--workers',type=int,default=3)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    named=set(a.names)
    jobs=[(j,json.loads((a.grid/(j[0]+'.json')).read_text())) for j in inputs(a.cache)
          if j[0] in named and not (a.out/(j[0]+'.json')).exists()]
    with mp.Pool(a.workers) as pool:
        for row in pool.imap_unordered(worker,jobs):
            (a.out/(row['event']+'.json')).write_text(json.dumps(row)+'\n')
