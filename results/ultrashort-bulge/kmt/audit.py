"""Stream, image-quality, timing and recurrence audit for *all* qualifying capped KMT grid fits."""
import argparse
import csv
import json
import multiprocessing as mp
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np

from freeze import HERE
from flux import stream_flux_scale
sys.path.insert(0,str(HERE.parent))
import fit


def read_streams(path, dia, year):
    # Evaluate DIAPL-only objects in the discovery season; products can span many years.
    from datetime import date
    start=(date(year,1,1)-date(1858,11,17)).days+.5
    with tarfile.open(path) as h:
        for m in h:
            if not m.isfile() or not m.name.endswith(('.diapl','.pysis')):continue
            cols=(0,1,2,4) if dia else (0,1,2,5)
            raw=h.extractfile(m).read()
            from io import BytesIO
            a=np.loadtxt(BytesIO(raw),usecols=cols,ndmin=2)
            cap=np.inf
            if not dia:
                if '_I.' in m.name:
                    b=np.loadtxt(BytesIO(raw),usecols=(1,3,4),ndmin=2)
                    base=stream_flux_scale(b)[0]
                    if np.isfinite(base):cap=1.15*base
                a[:,0]-=2400000;a[:,1]*=-1
            a=a[(a[:,0]>=start)&(a[:,0]<start+366)&np.isfinite(a[:,:3]).all(axis=1)&(a[:,2]>0)]
            if len(a)>=12:yield m.name,a,cap


def qualifying(r):
    out=[]
    if r['ps_te']<.5:out.append('ps')
    if r['fs_te']<.5 or r['fs_rho']>1:out.append('fs')
    return out


def shape(t,center,te,u,rho,engine):
    tau=(t-center)/te
    sep=np.hypot(u,tau)
    a=fit.lens.magnification(sep**2)-1.
    if rho:
        near=sep<max(2.,6*rho)
        if near.any():
            a[near]=np.asarray(engine.ESPLLightCurve([np.log(u),np.log(te),0.,np.log(rho)],(t[near]-center).tolist())[0])-1
    return a


def gain(x,y,e,mask,cap):
    if mask.sum()<3:return 0.
    x,y,w=x[mask],y[mask],e[mask]**-2
    x=x-np.average(x,weights=w)
    xx=np.dot(w,x*x)
    if xx<1e-15:return 0.
    xy=np.dot(w*x,y)
    fs=np.clip(xy/xx,0.,cap)
    return float(fs*(2*xy-fs*xx))


def diagnose(job):
    event,r,cat,path,dia=job
    try:
        engine=fit.engine(); year=int(event[4:8])
        datasets=list(read_streams(path,dia,year));out=[]
        for fam in qualifying(r):
            center=r[fam+'_t0']-2400000
            te=r[fam+'_te'];impact=r[fam+'_u0'];rho=r['fs_rho'] if fam=='fs' else 0.
            # FWHM of the template excess (a visible timescale for both PSPL and FSPL).
            from scipy.optimize import brentq
            def excess(tau):
                sep=np.hypot(impact,tau)
                return (engine.ESPLMag2(sep,rho) if rho else (sep**2+2)/(sep*np.sqrt(sep**2+4)))-1.
            half=te*brentq(lambda tau:excess(tau)-.5*excess(0.),0.,40*max(1.,impact,rho))
            totals=defaultdict(float);allg=qc=vm=0.;before=after=npeak=nhi=nbad=0;rec=set();sites=set()
            for stream,a,cap in datasets:
                t,f,e,seeing=a.T;quiet=abs(t-center)>max(5.,3*half)
                if quiet.sum()<10:quiet=np.ones(len(t),bool)
                floor=np.quantile(e[quiet],.05)
                err=np.maximum(e,floor);base=np.median(f[quiet])
                zz=(f[quiet]-base)/err[quiet];err*=max(1.,1.4826*np.median(abs(zz-np.median(zz))))
                # Include invalid and high-seeing exposures for first estimate.
                good=quiet&(seeing>0)&np.isfinite(seeing)
                threshold=np.quantile(seeing[good],.95) if good.sum()>=10 else np.inf
                q=(seeing>0)&(seeing<=threshold)
                a_model=shape(t,center,te,impact,rho,engine)
                near=abs(t-center)<=half
                # Compare to locally fitted constant, add a 10d baseline to define zero level.
                win=abs(t-center)<=max(10.,5*half)
                if '_V.' in stream:
                    # Ignore V-template support when all V points miss half maximum;
                    # the noise can spuriously favor a huge amplitude on a tiny template tail.
                    if near.any():vm+=gain(a_model,f,err,win&q,cap)
                    continue
                before+=int((near&(t<center)).sum());after+=int((near&(t>=center)).sum())
                npeak+=int(near.sum());z=(f-base)/err
                nhi+=int((near&(z>5)).sum());nbad+=int((near&(z>5)&~q).sum())
                allg+=gain(a_model,f,err,win,cap)
                good_gain=gain(a_model,f,err,win&q,cap);qc+=good_gain
                site=stream[:4]
                totals[site]+=good_gain
                if (near&q&(z>3)).sum()>=2:sites.add(site)
                # Recurrence flags use off-peak nights with at least three bright high-quality
                # samples, forming at least half the data in a one-day bin.
                mask=q&(abs(t-center)>max(12.,3*half))
                if mask.any():
                    days=np.floor(t[mask]).astype(int);first=int(days.min());days-=first
                    count=np.bincount(days);hits=np.bincount(days,weights=z[mask]>5,minlength=len(count))
                    for day in np.flatnonzero((hits>=3)&(2*hits>=count)):
                        rec.add((site,int(day+first)))
            other=sorted(totals.values(),reverse=True)
            out.append(dict(event=event,model=fam,te_d=te,rho=rho,t0=center,t0_minus_catalog_d=round(center-cat,5) if cat else None,
                            halfwidth_halfmax_hours=round(24*half,3),n_in_fwhm=npeak,n_before=before,n_after=after,
                            n_gt5=nhi,n_gt5_bad_seeing=nbad,sites_with_two_gt3='|'.join(sorted(sites)),
                            raw_gain=round(allg,3),quality_gain=round(qc,3),second_site_quality_gain=round(other[1] if len(other)>1 else 0.,3),
                            v_quality_gain=round(vm,3),offpeak_site_night_flags=len(rec),
                            site_quality_gains='|'.join(f'{s}:{v:.2f}' for s,v in sorted(totals.items())),status='ok'))
        return event,out
    except (ValueError,OSError,KeyError,TypeError,IndexError) as exc:
        return event,[dict(event=event,status='error',error=str(exc))]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--grid',type=Path,required=True);p.add_argument('--workers',type=int,default=2)
    args=p.parse_args();checkpoint=args.cache/'detector-audit';checkpoint.mkdir(exist_ok=True)
    ix={r['event']:r for r in map(json.loads,(args.cache/'index.jsonl').read_text().splitlines())}
    jobs=[]
    for pp in args.grid.glob('KMT*.json'):
        event=pp.stem
        if (checkpoint/pp.name).exists():continue
        r=json.loads(pp.read_text())
        if r.get('status')!='ok' or not qualifying(r):continue
        path=args.cache/(event+'.tar.gz');dia=not path.exists()
        if dia:path=args.cache/'fallback'/(event+'.tar.gz')
        if path.exists():jobs.append((event,r,ix[event]['t0'],path,dia))
    with mp.Pool(args.workers) as pool:
        for event,rows in pool.imap_unordered(diagnose,jobs):
            (checkpoint/(event+'.json')).write_text(json.dumps(rows)+'\n')
    rows=[r for p in sorted(checkpoint.glob('*.json')) for r in json.loads(p.read_text())]
    if rows:
        with (HERE/'short-lead-audit.tsv').open('w') as h:
            w=csv.DictWriter(h,list(dict.fromkeys(k for r in rows for k in r)),delimiter='\t',lineterminator='\n')
            w.writeheader();w.writerows(rows)


if __name__=='__main__':
    main()
