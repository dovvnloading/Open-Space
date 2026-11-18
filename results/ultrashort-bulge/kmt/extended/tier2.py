"""Multi-season and alias checks for five additional large-source flags."""
import argparse
import csv
import hashlib
import io
import json
import tarfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution,minimize

from sightlines import HERE,KMT,freeze_event
from recurrence import analyze,kernel,score,binned
from support import engine,read,median_abs_std
from pathlib import Path

GEOMETRY={
    'KMT-2024-BLG-1762':(60512.4790198,2.59635779,.00200005,15.),
    'KMT-2025-BLG-1797':(60945.5037771,13.233534,1.74959814,3.24936655),
    'KMT-2025-BLG-2489':(60942.7069485,16.917490,.002,2.61741915),
    'KMT-2024-BLG-1152':(60462.9155275,25.1612734,.002,1.66619544),
    'KMT-2024-BLG-3413':(60515.3203847,18.1183740,.002,2.46868182),
}


def fallback_overlap(cache,archives):
    event='KMT-2024-BLG-1762'
    file=archives/'fallback'/(event+'.tar.gz')
    rec=next(r for r in json.loads((KMT/'fallback-receipts.json').read_text()) if r['event']==event)
    if hashlib.sha256(file.read_bytes()).hexdigest()!=rec['curve']['sha256']:raise ValueError(file)
    src=pd.read_csv(cache/event/'BLG31M0105.058473.cache')
    details=[]
    key={'A':'SSO','C':'BLG','S':'SAO'}
    with tarfile.open(file) as archive:
        for m in archive:
            if not m.isfile() or not m.name.endswith('.diapl'):continue
            arr=np.loadtxt(archive.extractfile(m))
            df=pd.DataFrame(arr[:,:4],columns=['date','fallback_flux','fallback_error','fallback_chi2'])
            ok=(df.fallback_error>0)&(df.fallback_chi2>0)
            matched=df[ok].merge(src[src.site==key[m.name[3]]],on='date')
            details.append(dict(stream=m.name,valid_input=int(sum(ok)),matched_epochs=len(matched),
                                identical_flux=int(sum(np.isclose(matched.flux,matched.fallback_flux,atol=.001,rtol=0))),
                                by_year={str(y):len(g) for y,g in matched.groupby('year')}))
    return dict(receipt=rec['curve'],streams=details)


def duration_coverage(df,ge,year):
    t0,te,u,rho=ge
    # Intersections with the nominal source limb (rectilinear trajectory).
    half=te*np.sqrt(max(0.,rho**2-u**2))
    out={}
    for s,group in df[df.year==year].groupby('site'):
        k=group[(group.err>0)&(group.chi2>0)&(group.see>0)];t=k.date.to_numpy()
        out[s]=dict(n=len(t),first=round(t.min(),5),last=round(t.max(),5),
                    n_pre_ingress=int(sum(t<t0-half)),n_in_transit=int(sum((t>=t0-half)&(t<=t0+half))),n_post_egress=int(sum(t>t0+half)),
                    last_minus_predicted_egress=round(float(t.max()-t0-half),4))
    return out


def optical_amplitudes(archive,ge):
    center=ge[0];sc,tw,shape=kernel(ge)
    streams=[]
    with tarfile.open(archive) as tar:
        for member in tar:
            if not member.isfile() or not member.name.endswith('.pysis'):continue
            a=np.loadtxt(tar.extractfile(member),usecols=(0,1,2),ndmin=2)
            t,y,e=a[np.isfinite(a).all(axis=1)&(a[:,2]>0)].T
            t-=2400000+center;y=-y
            off=abs(t)>1.5*sc
            rec=dict(stream=member.name,n=len(t),n_within_rho_te=int(sum(abs(t)<sc)),n_outside_1p5_rho_te=int(sum(off)))
            if sum(off)>=5 and len(t)>=15:
                res=(y[off]-np.median(y[off]))/e[off]
                noise=max(1.,np.median(abs(res-np.median(res)))/.67449)
                w=(noise*e)**-2;h=np.interp(abs(t),tw,shape,right=0)
                h0=h-np.average(h,weights=w);den=np.sum(w*h0**2)
                if den>0:
                    amp=np.sum(w*h0*y)/den
                    residual=y-amp*h-np.average(y-amp*h,weights=w)
                    day=np.floor(t+center).astype(int);by=np.bincount(day-day.min(),weights=w*h0*residual)
                    sem=np.sqrt(np.sum(by**2))/den
                    if sem>0:rec.update(peak_excess_counts=round(float(amp),3),night_cluster_se=round(float(sem),3),
                                        S_N=round(float(amp/sem),3))
            streams.append(rec)
    return streams


def ogle_model(cache):
    cat=next(r for r in csv.DictReader((KMT.parent/'ogle/photometry.csv').open()) if r['event']=='2025-BLG-1203')
    raw,rec=read(cache,'ogle1203',cat['url'])
    if hashlib.sha256(raw).hexdigest()!=cat['sha256']:raise ValueError('OGLE raw hash drift')
    date,mag,em=np.loadtxt(io.BytesIO(raw),usecols=(0,1,2)).T;date-=2400000
    t0,te,u,rho=GEOMETRY['KMT-2025-BLG-1797']
    other=(date<60700)&(em>0)
    b=np.median(mag[other]);f=10**(-.4*(mag-b))
    err=np.log(10)*.4*f*em
    # Difference between individual observing-season baseline medians is retained in the report;
    # only 2025 exposures enter the fit. Calibrate uncertainties on earlier baseline photometry.
    z=(f[other]-np.median(f[other]))/err[other]
    scale=max(1.,median_abs_std(z))
    sel=(date>60700)&(date<61000)&(err>0)
    t,y,e=date[sel]-t0,f[sel],scale*err[sel];w=e**-2;v=engine()
    def flux(A):
        h=A-1.;hc=h-np.average(h,weights=w);den=np.sum(w*hc**2)
        amp=float(np.clip(np.sum(w*hc*y)/den,0.,1.15)) if den>1e-14 else 0.
        off=np.average(y-amp*h,weights=w)
        return float(np.sum(w*(y-off-amp*h)**2)),amp
    def chi(x,kind):
        center,logt,logu=x[:3];T,U=np.exp([logt,logu])
        if kind=='ps':
            s=((t-center)/T)**2+U**2;a=(s+2)/np.sqrt(s*(s+4))
        else:
            a=np.array(v.ESPLLightCurve([logu,logt,center,x[3]],t.tolist())[0])
        return flux(a)
    bounds=[(-45.,45.),(np.log(.1),np.log(300.)),(np.log(.002),np.log(5.))]
    models={}
    for kind,rr in [('ps',None),('fs',(np.log(.008),np.log(15.))),('fs_small',(np.log(.008),0.))]:
        bs=bounds+([rr] if rr else [])
        fun=lambda x:chi(x,kind)[0]
        seed=differential_evolution(fun,bs,seed=3641,popsize=8,maxiter=85,tol=.00005,polish=False)
        starts=[seed.x,[0.,np.log(te),np.log(u)]+([np.clip(np.log(rho),*rr)] if rr else [])]
        best=min([minimize(fun,x,method='Nelder-Mead',bounds=bs,options=dict(maxiter=1600,xatol=.0001,fatol=.001)) for x in starts],key=lambda x:x.fun)
        x=best.x;models[kind]=dict(chi2=round(float(best.fun),3),t0=round(float(t0+x[0]),4),
                       te=round(float(np.exp(x[1])),4),u0=round(float(np.exp(x[2])),4),
                       rho=round(float(np.exp(x[3])),4) if rr else 0.,F_source_over_baseline=round(chi(x,kind)[1],4))
    fixed=chi([0.,np.log(te),np.log(u),np.log(rho)],'fs')
    models['fixed_KMT_FS']=dict(chi2=round(fixed[0],3),F_source_over_baseline=round(fixed[1],4))
    half=te*np.sqrt(rho**2-u**2)
    return dict(receipt=rec,n_full=len(date),n_2025=len(t),baseline_I=float(b),baseline_error_scale=round(scale,4),
                n_in_nominal_crossing_2025=int(sum(abs(t)<half)),n_post_egress_2025=int(sum(t>half)),max_time_2025=float(date[sel].max()),
                constant_chi2=round(flux(np.ones_like(t))[0],3),models=models)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True);p.add_argument('--archives',type=Path,required=True)
    args=p.parse_args();args.cache.mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(2) as pool:
        sightlines=list(pool.map(lambda e:freeze_event(e,args.cache),GEOMETRY))
    (HERE/'tier2-sightlines.json').write_text(json.dumps(sightlines,indent=2)+'\n')
    out=[]
    for r in sightlines:
        e=r['event'];ge=GEOMETRY[e];yr=int(e.split('-')[1])
        d=analyze(e,r,args.cache,args.archives,geometry=ge,fallback=e.endswith('1762'))
        file=args.cache/e/(r['reference_stars'][0]['lcid']+'.cache');df=pd.read_csv(file)
        d['duration_coverage']=duration_coverage(df,ge,yr)
        if e.endswith(('1152','3413')):d['optical_amplitudes']=optical_amplitudes(args.archives/(e+'.tar.gz'),ge)
        out.append(d)
    (HERE/'tier2-results.json').write_text(json.dumps(dict(events=out,OGLE_1203=ogle_model(args.cache),
        fallback_overlap=fallback_overlap(args.cache,args.archives)),indent=2)+'\n')


if __name__=='__main__':main()
