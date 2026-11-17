"""Independent DIA, V-band and same-season diagnostics for KMT shortlist events."""
import argparse
import csv
import json
import sys
import tarfile
from pathlib import Path
from datetime import date
import numpy as np

from freeze import get, HERE
from controls import LEADS
from flux import baseline_caps
sys.path.insert(0,str(HERE.parent))
import fit


def stable_scale(z):
    return max(1.,1.4826*np.median(np.abs(z-np.median(z))))


def geometric_coverage(t, hours):
    # Centers yielding at least three distinct one-second timestamps in the centered window.
    t=np.unique(np.rint(t*86400)/86400)
    if len(t)<3 or t[-1]<=t[0]:
        return 0.
    h=hours/48
    a,b=np.maximum(t[2:]-h,t[0]),np.minimum(t[:-2]+h,t[-1])
    ok=b>a
    if not ok.any():
        return 0.
    a,b=a[ok],b[ok]
    end=np.r_[t[0],np.maximum.accumulate(b)[:-1]]
    return round(float(np.maximum(b-np.maximum(a,end),0.).sum()/(t[-1]-t[0])),5)


def measure(data, center, year, engine, model):
    t,flux,e=data.T
    good=np.isfinite(data).all(axis=1)&(e>0)&(t>50000)&(t<70000)
    t,flux,e=t[good],flux[good],e[good]
    if not len(t):return None
    # Separate season baselines, variable error scale, recurring positive nights.
    # Interseason gaps are much longer than duration; use civil year via day boundaries.
    offsets=np.array([(date(y,1,1)-date(1858,11,17)).days+.5 for y in range(2014,2029)])
    yr=2014+np.searchsorted(offsets,t,side='right')-1
    med=np.zeros(len(t));err=e.copy()
    for y in np.unique(yr):
        m=yr==y; quiet=m&(abs(t-center)>max(12.,3*min(20.,model['fs_te'])))
        if quiet.sum()<5: quiet=m
        med[m]=np.median(flux[quiet])
        err[m]=np.maximum(err[m],np.quantile(e[m],.05))
        err[m]*=stable_scale((flux[quiet]-med[quiet])/err[quiet])
    z=(flux-med)/err
    days=np.floor(t).astype(int)
    n=np.bincount(days)
    n5=np.bincount(days,weights=z>5,minlength=len(n))
    # Clusters are flags for review; seeing and weather can correlate an entire night.
    dates=np.flatnonzero((n>=3)&(n5>=3)&(2*n5>=n))
    offdays=dates[abs(dates+.5-center)>12]
    u=model['fs_u0'];te=model['fs_te'];rho=model['rho']
    tau=(t-center)/te
    sep=np.hypot(u,tau)
    magn=fit.lens.magnification(sep**2)-1
    close=sep<max(2.,6*rho)
    if close.any():
        magn[close]=np.asarray(engine.ESPLLightCurve([np.log(u),np.log(te),0.,np.log(rho)],(t[close]-center).tolist())[0])-1
    test=(yr==year)&(abs(t-center)<20.)
    near=abs(t-center)<=max(.15,te*max(1.,rho))
    def weighted_fit():
        if test.sum()<5:return (0.,None,None)
        x,y,w=magn[test],flux[test],err[test]**-2
        xbar=np.average(x,weights=w);ybar=np.average(y,weights=w)
        xx=np.dot(w,(x-xbar)**2)
        if xx<1e-15:return (0.,None,None)
        fs=np.dot(w*(x-xbar),y-ybar)/xx
        return float(max(fs,0.)**2*xx),float(fs),float(xx**-.5)
    gain,ampl,ampl_err=weighted_fit()
    return dict(n=len(t),first_t=round(t.min(),5),last_t=round(t.max(),5),observed_years='|'.join(map(str,np.unique(yr))),
                n_near=int(near.sum()),n_before=int((t<center-12).sum()),n_after=int((t>center+12).sum()),
                peak_gain=round(gain,2),source_flux=round(ampl,4) if ampl is not None else None,
                source_flux_sigma=round(ampl_err,4) if ampl_err is not None else None,
                n_recurring_nights=len(offdays),off_dates='|'.join(map(str,offdays)),peak_center=round(center,5))


def dump(path,rows):
    if not rows:return
    with path.open('w') as h:
        w=csv.DictWriter(h,fieldnames=list(dict.fromkeys(k for r in rows for k in r)),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--refined',type=Path,help='Override frozen stream-geometry.tsv with fit JSONs');p.add_argument('--event',action='append',default=[])
    args=p.parse_args()
    diaroot=args.cache/'diapl';diaroot.mkdir(exist_ok=True)
    leads=LEADS|set(args.event)
    frozen={}
    if (HERE/'stream-geometry.tsv').is_file():
        with (HERE/'stream-geometry.tsv').open() as h:
            frozen={r['event']:{k:float(v) for k,v in r.items() if k!='event'} for r in csv.DictReader(h,delimiter='\t')}
    leads|=set(frozen);geometries=[]
    ix={r['event']:r for r in map(json.loads,(args.cache/'index.jsonl').read_text().splitlines())}
    engine=fit.engine(); results=[];cadence=[]; receipts=[]
    for name in sorted(leads):
        refine=args.refined/(name+'.json') if args.refined else None
        model=json.loads(refine.read_text()) if refine and refine.exists() else frozen.get(name,{})
        if not model.get('rho'):continue
        geometries.append({'event':name,**{k:model[k] for k in ('fs_t0','fs_u0','fs_te','rho')}})
        year=int(name[4:8]);center=model['fs_t0']-2400000
        streams=[]
        kmt=args.cache/(name+'.tar.gz')
        if not kmt.exists():continue
        caps=baseline_caps(kmt)
        diafile=diaroot/(name+'.tar.gz')
        if not diafile.exists():
            try:
                data,r=get(ix[name]['url'].replace('pysis','diapl'))
                diafile.write_bytes(data);(diaroot/(name+'.json')).write_text(json.dumps(r))
            except OSError:pass
        if diafile.exists():receipts.append(dict(event=name,**json.loads((diaroot/(name+'.json')).read_text())))
        for archive in [kmt,diafile]:
            if not archive.exists():continue
            with tarfile.open(archive) as h:
                for m in h:
                    if not m.isfile() or not m.name.endswith(('.diapl','.pysis')):continue
                    a=np.loadtxt(h.extractfile(m),usecols=(0,1,2),ndmin=2)
                    if m.name.endswith('.pysis'):
                        a[:,0]-=2400000;a[:,1]*=-1
                        if '_I.' in m.name: streams.extend(a[(a[:,2]>0)&np.isfinite(a).all(axis=1),0])
                    row=measure(a,center,year,engine,model)
                    if row:
                        results.append(dict(event=name,stream=m.name,calibrated_I_baseline=caps.get(m.name),**row))
        t=np.sort(streams)
        cadence.append(dict(event=name,n=len(t),first_t=t.min(),last_t=t.max(),**{
            f'center_fraction_{h}h_atleast3':geometric_coverage(t,h) for h in (3,6,12,24)}))
    dump(HERE/'stream-geometry.tsv',geometries)
    dump(HERE/'stream-vetting.tsv',results)
    dump(HERE/'shortlist-cadence.tsv',cadence)
    dump(HERE/'dia-receipts.tsv',receipts)


if __name__=='__main__':
    main()
