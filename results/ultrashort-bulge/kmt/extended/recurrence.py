"""Profile fixed event shapes against public KMT field curves over multiple seasons."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import numpy as np
import pandas as pd

from support import HERE,KMT,engine,median_abs_std,cadence_colour
from sightlines import catalog

# Frozen public pySIS geometry (the fixed-shape amplitude test is independent).
GEOMETRY = {
    'KMT-2024-BLG-0824':(60441.5737763,7.8999843,.34659064,2.0610926),
    'KMT-2024-BLG-3278':(60422.6748584,8.15917179,.00200153,2.58195269),
    'KMT-2025-BLG-1495':(60930.3743737,33.8324713,.18679515,2.898967),
    'KMT-2025-BLG-0220':(60884.8592940,.10855826,.002,5.0515451),
}


def binned(df,ge):
    mid,te,_,rho=ge
    cross=te*rho
    # One observing night for long peaks and 0.05 d for the hours-long alert.
    interval=1. if te>1 else .05
    df=df[np.isfinite(df[['date','flux','err','chi2','see','year']]).all(axis=1)&(df.err>0)&(df.chi2>0)&(df.see>0)].copy()
    df['bin']=np.floor(df.date/interval).astype(int)
    obs={}
    for (site,yr),group in df.groupby(['site','year']):
        off=group[abs(group.date-mid)>1.5*cross]
        if len(off)<12:continue
        group=group[group.see<=np.quantile(off.see,.95)]
        n=group.groupby('bin').agg(t=('date','median'),y=('flux','median'),e=('err','median'),n=('date','size'))
        if te>1:n=n[n.n>=2]
        ctrl=n[abs(n.t-mid)>1.5*cross]
        if len(ctrl)<8:continue
        noise=max(median_abs_std(ctrl.y.to_numpy()),1e-6)
        e=np.hypot(noise,n.e.to_numpy()/np.sqrt(n.n.to_numpy()))
        # Isolated extreme nights can dominate large historical deviations. Cap flux to
        # ±8 season/site scatters on every date, including the actual event.
        y=np.clip(n.y.to_numpy()-np.median(ctrl.y),-8*noise,8*noise)
        obs.setdefault(int(yr),{})[site]=(n.t.to_numpy(),y,e)
    return obs


def score(obs,c,scale,tw,sh):
    out={}
    for site,(t,y,e) in obs.items():
        # Sorted bins permit reasonably fast coverage-window searches.
        lo,hi=np.searchsorted(t,[c-tw[-1],c+tw[-1]])
        dt,y0,e0=abs(t[lo:hi]-c),y[lo:hi],e[lo:hi]
        if len(dt)==0 or sum(dt<scale)<2 or sum(dt>1.35*scale)<4:
            continue
        a=np.interp(dt,tw,sh)
        w=e0**-2
        center=a-np.average(a,weights=w)
        den=np.sum(w*center**2)
        if den<1.e-12:continue
        amp=np.sum(w*center*y0)/den
        out[site]=[round(amp,2),round(amp*np.sqrt(den),2),int(sum(dt<scale)),int(sum(dt>1.35*scale))]
    # The second site's signed S/N suppresses huge single-site systematics.
    scores=sorted([v[1] for v in out.values()],reverse=True)
    return round(scores[1],2) if len(scores)>=2 else None,out


def kernel(ge):
    _,te,u,rho=ge;scale=te*rho
    tw=np.linspace(0,3.5*scale,1536)
    l=engine();sh=np.array([l.ESPLMag2(np.hypot(x/te,u),rho)-1 for x in tw])
    return scale,tw,sh/sh[0]


def analyze(event,records,cache,archives,geometry=None,fallback=False):
    ge=GEOMETRY[event] if geometry is None else geometry;c,te,u,rho=ge
    yr=int(event.split('-')[1]);scale,tw,sh=kernel(ge)
    all_obs=[]
    for r in records['reference_stars']:
        f=cache/event/(r['lcid']+'.cache')
        if hashlib.sha256(f.read_bytes()).hexdigest()!=r['receipt']['sha256']:raise ValueError(f)
        all_obs.append(binned(pd.read_csv(f),ge))
    centers=[('geometry',yr,c)]
    results=[]
    obs=all_obs[0]
    obs_stats={}
    for year,sites in obs.items():
        obs_stats[year]={s:len(d[0]) for s,d in sites.items()}
        first=min(v[0].min() for v in sites.values());last=max(v[0].max() for v in sites.values())
        # Grid spacing resolves the shortest transit and avoids oversized searches for long pulses.
        scans=np.arange(np.floor(first),np.ceil(last),min(4.,scale/4))
        best=None;eligible=0
        for x in scans:
            if abs(x-c)<=2*scale:continue
            s,_=score(sites,x,scale,tw,sh)
            if s is None:continue
            eligible+=1
            if best is None or s>best[0]:best=s,x
        if best:
            centers.append((f'annual_max_{eligible}_centers',year,best[1]))
    row=catalog(event);alert=float(row['t_0']) if row['t_0']!='-' else None
    if alert is not None and abs(c-alert)>.5*scale:centers.append(('catalog_alert',yr,alert))
    for tag,year,x in centers:
        values=[]
        for r,ob in zip(records['reference_stars'],all_obs):
            sn,by=score(ob.get(year,{}),x,scale,tw,sh)
            values.append(dict(id=r['lcid'],distance_arcsec=round(r['dist'],3),second_site_sn=sn,sites=by))
        results.append(dict(kind=tag,year=year,center=round(float(x),5),position_fits=values))
    geom=dict(t0=c+2400000,te_or_rise_d=te,u0_or_fall_d=u,rho=rho)
    if fallback:
        colors=[]
    else:
        file=archives/(event+'.tar.gz')
        with (KMT/f'pySIS-{yr}.tsv').open() as f:
            r=next(r for r in csv.DictReader(f,delimiter='\t') if r['event']==event)
        if hashlib.sha256(file.read_bytes()).hexdigest()!=r['sha256']:raise ValueError(file)
        colors=cadence_colour(file,geom)
    return dict(event=event,fixed_geometry=geom,pySIS_color=colors,
                season_site_bin_counts=obs_stats,tests=results)


def ogle_1000(file):
    with (KMT.parent/'ogle/photometry.csv').open() as f:
        receipt=next(r for r in csv.DictReader(f) if r['event']=='2025-BLG-1000')
    if hashlib.sha256(file.read_bytes()).hexdigest()!=receipt['sha256']:raise ValueError(file)
    t,m,em=np.loadtxt(file,usecols=[0,1,2]).T
    ok=np.isfinite(t)&np.isfinite(m)&np.isfinite(em)&(em>0)
    t,m,em=t[ok]-2400000,m[ok],em[ok]
    med=np.median(m);flux=10**(-.4*(m-med));e=np.log(10)*.4*flux*em
    mid=GEOMETRY['KMT-2025-BLG-1495'][0]
    scale,tw,sh=kernel(GEOMETRY['KMT-2025-BLG-1495'])
    # Intercepts independent for the event season; remove the distinct narrow catalog alert.
    sel=(t>60700)&(t<61000)&(abs(t-60850.779)>5)
    h=np.interp(abs(t[sel]-mid),tw,sh)
    w=e[sel]**-2;ch=h-np.average(h,weights=w);den=np.sum(w*ch**2)
    amp=np.sum(w*ch*flux[sel])/den
    eamp=den**-.5
    # Cluster standard error by integer night using the same weighted regression.
    residual=flux[sel]-np.average(flux[sel]-amp*h,weights=w)-amp*h
    day=np.floor(t[sel]).astype(int)
    byday=np.bincount(day-day.min(),weights=w*ch*residual)
    clustered=np.sqrt((byday**2).sum())/den
    checks=[]
    for cen in (60850.779,mid):
        for d in [.2,1.,5.,30.,70.]:
            k=abs(t-cen)<d
            checks.append(dict(center=round(cen,5),radius_d=d,n=int(sum(k)),median_I=round(float(np.median(m[k])),4),
                               median_err=round(float(np.median(em[k])),4)))
    return dict(receipt=receipt,n=int(len(t)),baseline_I=float(med),fixed_long_shape=dict(n=int(sum(sel)),min_t=t[sel].min(),
                max_t=t[sel].max(),peak_fraction=amp,error_independent=eamp,error_night_clustered=clustered),
                local=checks)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archives',type=Path,required=True)
    p.add_argument('--cache',type=Path,required=True);p.add_argument('--ogle-1000',type=Path,required=True)
    args=p.parse_args();receipt=json.loads((HERE/'sightline-receipts.json').read_text())
    result=[analyze(r['event'],r,args.cache,args.archives) for r in receipt]
    (HERE/'recurrence.json').write_text(json.dumps(dict(events=result,OGLE_1000=ogle_1000(args.ogle_1000)),indent=2)+'\n')


if __name__=='__main__':main()
