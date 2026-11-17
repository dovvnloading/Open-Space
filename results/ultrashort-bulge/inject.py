"""Paired injections into off-alert windows at real epochs. Calendar-uniform peak trials.

These trials condition on existing alert sightlines, a supplied peak epoch and an
imposed peak signal-to-noise per stream. They cannot calibrate survey discovery rates.
"""
import argparse
import csv
import hashlib
from pathlib import Path
import numpy as np
import fit
from run_grid import inputs,read_curve


def key(name):
    return int.from_bytes(hashlib.sha256(('ultrashort-v1 '+name).encode()).digest()[:8],'little')


def observations(job, centers):
    t,y,e,g,diff=read_curve(job);t,y,e,g=fit.clean(t,y,e,g)
    if diff:e=fit.reported_error_floor(e,g)
    if not len(t):return
    # Calibrate off-event noise separately per stream on the survey's own scale.
    level=np.zeros(max(g)+1);noise=level.copy()
    for j in range(len(level)):
        ok=(g==j)&(abs(t-job[2])>max(3.,min(60.,2*abs(job[3])))) if np.isfinite(job[2]+job[3]) else (g==j)
        if not np.any(ok):ok=g==j
        level[j]=np.median(y[ok]);noise[j]=np.median(e[ok])
    for center in centers:
        use=abs(t-center)<4
        yield center,t[use],y[use],e[use],g[use],level,noise,diff


def run(job,number):
    name=job[0];yr=2000+int(name[4:6]) if name.startswith('Gaia') else int(name.split('-')[1])
    start=2451544.5+int((np.datetime64(f'{yr}-03-01')-np.datetime64('2000-01-01')).astype(int))
    rng=np.random.default_rng(key(name));centers=[]
    # Exclude the known event's neighborhood identically for all trial widths.
    while len(centers)<number:
        time=start+rng.uniform(0.,245.)
        if not np.isfinite(job[2]+job[3]) or abs(time-job[2])>max(3.,min(60.,2*abs(job[3]))):
            centers.append(time)
    vb=fit.engine()
    for center,t,y,e,g,base,noise,diff in observations(job,centers):
        for te in (.03,.1,.3):
            for rho in (0.,2.):
                w=te*max(1.,rho)
                before=np.unique(t[(t<center)&(t>=center-w)])
                after=np.unique(t[(t>center)&(t<=center+w)])
                both=min(len(before),len(after))>=2
                r=dict(event=name,t0=round(center,6),te_d=te,rho=rho,peak_snr=12,two_each_side=int(both),
                       paired_pass=0,control_pass=0)
                if not both:
                    yield r;continue
                groups=np.unique(g)
                # One weighted baseline constraint per participating stream.
                sd=noise[groups]/np.sqrt(25.)
                ta=np.r_[t,np.full(len(groups),center+150.)]
                ya=np.r_[y,base[groups]];ea=np.r_[e,sd];ga=np.r_[g,groups]
                u0=.35
                u=np.hypot((t-center)/te,u0)
                if rho:
                    excess=np.array([vb.ESPLMag2(v,rho)-1 for v in u]);peak=vb.ESPLMag2(u0,rho)-1
                else:
                    excess=fit.lens.magnification(u*u)-1;peak=fit.lens.magnification(u0*u0)-1
                y_inj=ya.copy();y_inj[:len(t)]+=12*noise[g]*excess/peak
                def measurement(yy):
                    s=fit.screen(ta,yy,ea,ga,center,np.nan,yr,diff)
                    family=max(('fs','ps'),key=lambda k:s[k+'_gain'])
                    ok=(abs(s[family+'_t0']-center)<=w and (s[family+'_te']<.5 or (family=='fs' and s['fs_rho']>1))
                        and s[family+'_gain']>=s['bin_gain'])
                    return s[family+'_gain'],ok
                control,ctrl=measurement(ya);injected,success=measurement(y_inj)
                r.update(control_gain=round(control,2),injected_gain=round(injected,2),
                         paired_pass=int(success and injected-control>=80),control_pass=int(ctrl and control>=80))
                yield r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,default=Path.cwd())
    p.add_argument('--per-survey',type=int,default=24);p.add_argument('--seasons',type=int,default=8,
        help='trial centers per selected light curve');p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();fit.setup();all_jobs=list(inputs(a.cache))
    order=lambda j:key(j[0])
    with a.out.open('w') as out:
        writer=csv.DictWriter(out,['event','t0','te_d','rho','peak_snr','two_each_side','paired_pass',
                                   'control_pass','control_gain','injected_gain']);writer.writeheader()
        for survey in ('OGLE','KMT','MOA','PRIME','Gaia'):
            jobs=sorted([j for j in all_jobs if j[0].startswith(survey)],key=order)[:a.per_survey]
            for j in jobs:
                writer.writerows(run(j,a.seasons))
