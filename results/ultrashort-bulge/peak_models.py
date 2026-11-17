"""Profile competing peak shapes using the observed site/field streams.

Consistent local window and baseline summary for all families. Source flux can
be bounded by a multiple of the photometric baseline (KMT, OGLE), or unbounded
for difference-only data. Report the change under the chosen calibration.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import fit
from run_grid import inputs,read_curve
from kmt.flux import baseline_caps


def solve(job, seed, factor, binary_trials=12):
    name,path,t0,te=job
    vb=fit.engine(); vb.Tol=.001
    t,y,err,g,diff=read_curve(job)
    t,y,err,g=fit.clean(t,y,err,g)
    ref=seed['fs_t0']
    t=t-ref
    if diff: err=fit.reported_error_floor(err,g)
    caps=np.full(max(g)+1,np.inf)
    if name.startswith('KMT') and np.isfinite(factor):
        raw=np.asarray(list(baseline_caps(path).values()))
        if len(raw)!=len(caps):
            raise ValueError('stream alignment failed')
        caps=np.where(np.isfinite(raw),factor*raw,np.inf)
    if not diff:
        caps=np.array([factor*np.median(y[g==i]) for i in range(len(caps))])
    # A finite comparison window, with the out-of-window mean per site as a
    # single baseline measurement at A=1. Freeze clipping outside peak only.
    tt=[]; yy=[]; ee=[]; gg=[]
    for i in range(len(caps)):
        m=g==i
        out=m&(abs(t)>12)
        if out.sum()<4: out=m&(abs(t)>4)
        if out.sum()==0: continue
        z=(y[out]-np.median(y[out]))/err[out]
        scale=max(1.,np.median(abs(z))/.67448975)
        ok=abs(z)<6*scale
        wo=1/(scale*err[out][ok])**2
        if not len(wo): continue
        near=m&(abs(t)<=4)
        tt.extend(t[near]);yy.extend(y[near]);ee.extend(err[near]*scale);gg.extend([i]*near.sum())
        tt.append(1e6);yy.append(np.average(y[out][ok],weights=wo));ee.append(1/np.sqrt(wo.sum()));gg.append(i)
    t,y,e,g=map(np.asarray,(tt,yy,ee,gg))
    peak=t<1e5
    wt=1/e**2
    wg=np.bincount(g,weights=wt,minlength=len(caps))
    ym=np.bincount(g,weights=wt*y,minlength=len(caps))/np.maximum(wg,1e-100)
    fy=y-ym[g]
    def score(v,kind,ret=False):
        dt,lte,lu=v[:3];te,u=np.exp([lte,lu])
        amp=np.zeros(len(t));tp=t[peak]
        if kind=='ps': amp[peak]=fit.lens.magnification(((tp-dt)/te)**2+u*u)-1
        elif kind=='fs':
            amp[peak]=np.array(vb.ESPLLightCurve([lu,lte,dt,v[3]],tp.tolist())[0])-1
        elif kind=='bin':
            # VB's binary trajectories have their own axis sign convention;
            # both orientations are searched from independent angle starts.
            amp[peak]=np.array(vb.BinaryLightCurve([v[4],v[5],u,v[6],v[3],lte,dt],tp.tolist())[0])-1
        elif kind=='exp':
            amp[peak]=np.exp(np.clip(np.where(tp<dt,(tp-dt)/te,-(tp-dt)/u),-745,0))
        if not np.isfinite(amp).all(): return 1e20
        ax=amp-(np.bincount(g,weights=wt*amp,minlength=len(caps))/np.maximum(wg,1e-100))[g]
        xx=np.bincount(g,weights=wt*ax**2,minlength=len(caps))
        fs=np.clip(np.bincount(g,weights=wt*ax*fy,minlength=len(caps))/np.maximum(xx,1e-200),0.,np.inf if kind=='exp' else caps)
        chi=np.sum(wt*(fy-fs[g]*ax)**2)
        if ret:return fs.tolist()
        return chi
    rng=np.random.default_rng(623)
    base=[(-2,2),(np.log(.01),np.log(150)),(np.log(.002),np.log(12))]
    results={}
    for kind in ['ps','fs','exp','bin']:
        bounds=base[:]
        if kind in ('fs','bin'): bounds+=[(np.log(.008),np.log(15))]
        if kind=='bin':bounds += [(np.log(.3),np.log(3.5)),(np.log(.0003),0),(-7,7)]
        initial=[]
        for family in ('ps','fs'):
            start=[seed[family+'_t0']-ref,np.log(seed[family+'_te']),np.log(seed[family+'_u0'])]
            if kind in ('fs','bin'): start += [np.log(seed['fs_rho'])]
            initial.append(start)
        if kind=='exp': initial += [[0,-2,-1],[0,-1,-2]]
        if kind=='fs':
            initial += [list(results['ps']['vector'])+[np.log(.15)],[0,-.6,-1.,-.5],[0,-1,0,0]]
        if kind=='bin':
            initial=[]
            for trial in range(binary_trials):
                s,q,ang=(seed['bin_s'],seed['bin_q'],seed['bin_angle']) if trial==0 else (np.exp(rng.uniform(-.65,.65)),np.exp(rng.uniform(-5,0)),rng.uniform(-np.pi,np.pi))
                initial.append([rng.normal(0,.15), np.log(max(.3,min(seed['bin_te'],8.)))+rng.normal(0,.7),
                                np.log(rng.uniform(.07,.8)),np.log(.035),np.log(s),np.log(q),ang])
        best=None
        for v in initial:
            v=np.clip(v,*np.array(bounds).T)
            simp=np.tile(v,(len(v)+1,1))
            for k in range(len(v)):
                step=.06 if k==0 else .16
                simp[k+1,k]+=step if v[k]+step<bounds[k][1] else -step
            opt=minimize(score,v,args=(kind,),method='Nelder-Mead',bounds=bounds,
                         options={'initial_simplex':simp,'maxiter':420,'xatol':.0003,'fatol':.01})
            if best is None or opt.fun<best.fun: best=opt
        results[kind]={'chi2':round(best.fun,3),'t0':ref+best.x[0],'te_or_rise_d':np.exp(best.x[1]),
                       'u0_or_fall_d':np.exp(best.x[2]),'vector':best.x.tolist(),
                       'fs':score(best.x,kind,True),'converged':bool(best.success)}
        if kind in ('fs','bin'):results[kind]['rho']=np.exp(best.x[3])
        if kind=='bin':results[kind].update(s=np.exp(best.x[4]),q=np.exp(best.x[5]),alpha=best.x[6])
    return {'event':name,'n_peak_window':int(peak.sum()),'baseline_pseudopoints':int((~peak).sum()),
            'cap_factor':factor if np.isfinite(factor) else None,'caps':[d if np.isfinite(d) else None for d in caps],
            'reference_jd':ref,'models':results}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('event');p.add_argument('--cache',type=Path,default=Path.cwd())
    p.add_argument('--grid',type=Path,default=Path('staging/screen'));p.add_argument('--factor',type=float,default=1.15)
    p.add_argument('--trials',type=int,default=12);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();j=next(j for j in inputs(a.cache) if j[0]==a.event)
    seed=json.loads((a.grid/(a.event+'.json')).read_text())
    res=solve(j,seed,a.factor,a.trials);a.out.write_text(json.dumps(res,indent=2,allow_nan=False)+'\n')
