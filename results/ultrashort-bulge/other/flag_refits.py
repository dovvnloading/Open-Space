"""Continuously refit every non-KMT broad flag on robustly normalized peak windows.

Point lens, finite uniform source and two-sided exponential each fit a free
baseline and positive flux. The FS shape uses VBMicrolensing's interpolated ESPL
table. Local optimizations and sparse-data degeneracies limit the fitted ranking.
"""
import argparse
import csv
import gzip
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from run_grid import inputs,read_curve
import fit

WORK=Path('staging/flag-refits')
SEED=Path('staging/screen')
OUT=HERE/'all-flag-refits.csv.gz'


def window(job,s,calibration=None,cap_factor=1.15):
    event,path,cat_t,cat_te=job
    t,y,e,g,diff=read_curve(job)
    t,y,e,g=fit.clean(t,y,e,g)
    if diff:e=fit.reported_error_floor(e,g)
    origin=s['fs_t0']
    t=t-origin
    half=np.clip(4*max(s['ps_te'],s['fs_te']*max(1.,s['fs_rho'])),3.,130.)
    # Preserve each tested local peak when the PS and FS grid centers differ.
    dps=s['ps_t0']-origin
    left=min(0,dps)-half;right=max(0,dps)+half
    near=(t>=left)&(t<=right)
    if near.sum()<3:
        nearest=np.sort(np.minimum(abs(t),abs(t-dps)))[:min(6,len(t))]
        half=max(half,float(nearest[-1])+.01)
        left=min(0,dps)-half;right=max(0,dps)+half
        near=(t>=left)&(t<=right)
    far=(t<left-half/2)|(t>right+half/2)
    if far.sum()<3:far=~near
    native_far=far.copy()
    if far.sum()<3:
        # Single-season curves may show no secure baseline. Fit all native
        # fluxes; do not create an artificial zero-magnification constraint.
        e_scale=1.;base=np.quantile(y,.2);err_base=None
    else:
        mid=np.median(y[far])
        z=(y[far]-mid)/e[far]
        e_scale=max(1.,np.median(abs(z))/.67448975)
        ok=abs(z)<5*e_scale
        wf=1/(e[far][ok]*e_scale)**2
        base=np.average(y[far][ok],weights=wf)
        err_base=1/np.sqrt(sum(wf))
    e*=e_scale
    # Binning set by the narrowest competing grid timescale. A minimum
    # 6-minute bin retains hour-scale structure; large events allow up to 4.8 h.
    bin_d=np.clip(min(s['ps_te'],s['fs_te'])/14,.0042,.2)
    tn=t[near];yn=y[near];en=e[near]
    if len(tn)<3:raise ValueError('fewer than three observations inside both grid peak windows')
    group=np.unique(np.floor(tn/bin_d),return_inverse=True)[1]
    w=1/en**2
    totals=np.bincount(group,weights=w);n=np.bincount(group)
    tb=np.bincount(group,weights=w*tn)/totals
    yb=np.bincount(group,weights=w*yn)/totals
    # Inflate uncertainty in a bin with measurable extra dispersion rather
    # than allowing repeated inconsistent measurements to act independently.
    chi=np.bincount(group,weights=w*(yn-yb[group])**2)
    eb=np.sqrt(np.maximum(1.,chi/np.maximum(n-1,1))/totals)
    peak=np.ones(len(tb),bool)
    if err_base is not None:
        tb=np.r_[tb,1e6];yb=np.r_[yb,base];eb=np.r_[eb,err_base];peak=np.r_[peak,False]
    cap=np.inf if diff else max(1e-9,cap_factor*np.median(y))
    quiet_median=float(np.median(y[far])) if far.any() else float(np.median(y))
    if calibration is not None:
        fbase=float(calibration['additive_reference_counts'])+quiet_median
        if fbase>0:cap=cap_factor*fbase
    
    return (tb,yb,eb,peak,cap),(origin,half,bin_d,int(n.sum()),int(far.sum()),e_scale,quiet_median),(t,y,e,base)


def fit_models(job,seed,calibration=None,cap_factor=1.15):
    (t,y,e,pk,cap),meta,raw=window(job,seed,calibration,cap_factor)
    origin,hw,bins,n,off,escale,quiet_median=meta
    event=job[0];tn=t[pk]
    vb=fit.engine();vb.Tol=.0003
    wt=1/e**2;w=wt/wt.sum();ym=np.dot(w,y);std=y-ym
    null=np.dot(wt,std**2)
    def score(x,kind,ret=False):
        shift,lg,lu=x[:3];te,u=np.exp([lg,lu]);z=(tn-shift)/te
        if kind=='exp':
            signal=np.exp(np.clip(np.where(tn<shift,(tn-shift)/te,-(tn-shift)/u),-740,0))
        else:
            sep=np.hypot(u,z)
            signal=fit.lens.magnification(sep**2)-1
            if kind=='fs':
                r=np.exp(x[3])
                near=sep<max(6*r,2)
                if near.any():
                    signal[near]=np.array(vb.ESPLLightCurve([lu,lg,shift,x[3]],tn[near].tolist())[0])-1
        a=np.zeros(len(t));a[pk]=signal
        mean=np.dot(w,a);ax=a-mean
        flux=np.clip(np.dot(wt*ax,std)/max(np.dot(wt,ax**2),1e-200),0,np.inf if kind=='exp' else cap)
        r=(std-flux*ax)/e
        if ret:return r,flux,ym-flux*mean
        return np.dot(r,r)
    a=seed['ps_t0']-origin
    tlo=min(a,0)-max(.3,.43*hw);thi=max(a,0)+max(.3,.43*hw)
    bounds=[(tlo,thi),(np.log(.012),np.log(350)),(np.log(.008),np.log(10))]
    ps0=np.array([a,np.log(seed['ps_te']),np.log(seed['ps_u0'])])
    fs0=np.array([0,np.log(seed['fs_te']),np.log(seed['fs_u0'])])
    models={}
    def optimize(kind,starts,limits,budget):
        solutions=[]
        for ini in starts:
            ini=np.clip(ini,*np.array(limits).T)
            simplex=np.tile(ini,(len(ini)+1,1))
            for j in range(len(ini)):
                step=min(.6,max(.02,.08*min(seed['ps_te'],seed['fs_te']))) if j==0 else .22
                simplex[j+1,j]+=step if ini[j]+step<limits[j][1] else -step
            opt=minimize(score,ini,args=(kind,),method='Nelder-Mead',bounds=limits,
                         options=dict(initial_simplex=simplex,maxiter=budget,xatol=.0012,fatol=.018))
            solutions.append(opt)
        return min(solutions,key=lambda v:v.fun)
    models['ps']=optimize('ps',[ps0,fs0],bounds,115)
    models['fs']=optimize('fs',[np.r_[fs0,np.log(seed['fs_rho'])],np.r_[models['ps'].x,np.log(.18)],
                                np.r_[ps0,np.log(1.8)]],bounds+[(np.log(.03),np.log(15))],155)
    ew=max(.02,seed['fs_te']*max(1,seed['fs_rho']))
    models['exp']=optimize('exp',[[0,np.log(.22*ew),np.log(ew)],[0,np.log(ew),np.log(.22*ew)]],
                           bounds[:2]+[(np.log(.012),np.log(350))],110)
    out=dict(event=event,peak_window_halfwidth_d=round(hw,4),bin_d=round(bins,5),peak_epochs=n,peak_bins=int(pk.sum()),
             baseline_epochs=off,scatter_inflation=round(escale,3),constant_chi2=round(null,3),status='ok')
    if calibration is not None:
        fr=float(calibration['additive_reference_counts'])
        out.update(additive_reference_counts=fr,quiet_median_difference_counts=round(quiet_median,3),
                   baseline_counts=round(fr+quiet_median,3),cap_factor=cap_factor,
                   catalog_band=calibration['band'],zero_point=calibration['zero_point'],
                   implied_baseline_mag=round(float(calibration['zero_point'])-2.5*np.log10(fr+quiet_median),4) if fr+quiet_median>0 else '')
    for kind,sol in models.items():
        x=sol.x
        rr,flux,base=score(x,kind,True)
        hub=np.sum(np.where(abs(rr)<3,rr**2,6*abs(rr)-9))
        out.update({kind+'_chi2':round(sol.fun,3),kind+'_huber':round(hub,3),kind+'_t0':round(float(origin+x[0]),6),
                    kind+'_te_or_rise':round(float(np.exp(x[1])),6),kind+'_u0_or_fall':round(float(np.exp(x[2])),6),
                    kind+'_converged':bool(sol.success)})
        if kind=='fs':
            out['rho']=round(float(np.exp(x[3])),5);out['source_fraction']=round(float(cap_factor*flux/cap),5) if np.isfinite(cap) else ''
            tr,yr,er,baseline=raw
            r_out=(abs(tr-x[0])>max(10.,hw*1.25))
            # For isolated-peak ranking use a robust global noise floor and
            # require adjacent high points within one day to expose recurrence.
            noise=1.4826*np.median(abs(yr[r_out]-np.median(yr[r_out]))) if r_out.sum()>8 else 0.
            z=(yr-baseline)/np.hypot(er,noise)
            tp=np.exp(x[1])*max(1.,np.exp(x[3]),np.exp(x[2]));near=abs(tr-x[0])<min(hw,max(.12,1.5*tp))
            out.update(n_peak_before_3sigma=int(np.count_nonzero(near&(z>3)&(tr<x[0]))),
                       n_peak_after_3sigma=int(np.count_nonzero(near&(z>3)&(tr>=x[0]))))
            triple=np.minimum.reduce([z[:-2],z[1:-1],z[2:]])
            valid=r_out[:-2]&r_out[2:]&(tr[2:]-tr[:-2]<1.)
            ix=np.flatnonzero(valid)
            if len(ix):
                im=ix[np.argmax(triple[ix])]
                out['offpeak_triple_min_z']=round(float(triple[im]),2)
                out['offpeak_triple_hjd']=round(float(tr[im+1]+origin),5)
            out['catalog_peak_offset_d']=round(float(origin+x[0]-job[2]),4) if np.isfinite(job[2]) else ''
            out['peak_window_chi2_per_bin']=round(float(np.sum(rr[pk]**2)/max(pk.sum()-5,1)),3)
    out['fs_over_both']=round(min(out['ps_chi2'],out['exp_chi2'])-out['fs_chi2'],3)
    out['fs_over_both_huber']=round(min(out['ps_huber'],out['exp_huber'])-out['fs_huber'],3)
    return out


def fit_event(job):
    name=job[0]
    path=WORK/(name+'.json')
    if path.exists():return
    try:
        seed=json.loads((SEED/(name+'.json')).read_text())
        row=fit_models(job,seed)
    except Exception as err:
        row=dict(event=name,status='failed',error=repr(err))
    tmp=path.with_suffix('.tmp')
    tmp.write_text(json.dumps(row,allow_nan=False))
    tmp.replace(path)


def selected():
    events=[]
    for path in SEED.glob('*.json'):
        if path.name.startswith('KMT'):continue
        row=json.loads(path.read_text())
        if row.get('status')=='ok' and (row['ps_te']<.5 or row['fs_te']<.5 or row['fs_rho']>1):
            # Fit strong events and all 2026 zero-lag candidates first.
            rank=(min(row['ps_te'],row['fs_te'])>=.5,'2026' not in path.stem,
                  -min(row['fs_gain']-row['ps_gain'],row['fs_gain']-row['bin_gain']))
            events.append((rank,path.stem))
    return [n for _,n in sorted(events)]


def publish(names):
    from consolidate import write_gz
    result=[json.loads((WORK/(n+'.json')).read_text()) if (WORK/(n+'.json')).exists()
            else dict(event=n,status='missing checkpoint') for n in sorted(names)]
    write_gz(OUT,result)
    from collections import Counter
    print(len(result),Counter(r['status'] for r in result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--workers',type=int,default=4);p.add_argument('--publish-only',action='store_true')
    p.add_argument('--limit',type=int);args=p.parse_args()
    WORK.mkdir(parents=True,exist_ok=True)
    names=selected()
    if args.limit:names=names[:args.limit]
    if not args.publish_only:
        by_id={j[0]:j for j in inputs(Path.cwd()) if not j[0].startswith('KMT')}
        with ProcessPoolExecutor(args.workers) as workers:
            list(workers.map(fit_event,[by_id[n] for n in names if n in by_id],chunksize=1))
    publish(names)
