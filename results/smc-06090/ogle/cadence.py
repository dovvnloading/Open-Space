"""Validate the measured separation change at the real early and late observation times."""
import argparse
import json
from pathlib import Path
import numpy as np
from timing import OBJECT, derived, fit, obtain, profile, signed_phase

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',default='timing-inputs',type=Path)
    p.add_argument('--output',type=Path,default='cadence.json')
    p.add_argument('--draws',type=int,default=800)
    a=p.parse_args(); obtain(a.cache)
    row=next(l.split() for l in (a.cache/'ecl.dat').read_text().splitlines() if l.split()[0]==OBJECT)
    period,epoch=map(float,row[3:5])
    t,m,e=np.loadtxt(a.cache/'ogle4.dat').T
    rng=np.random.default_rng(803122)
    x=[]
    for selected in [t<2457000,t>=2459700]:
        ok=selected & np.isfinite(t+m+e) & (e>0)&(e<=.10)
        ph=signed_phase((t[ok]-epoch)/period); mag,err=m[ok],e[ok]
        best=fit(ph,mag,err,[np.percentile(mag,25),.4,0,.023,.4,.51,.023], 'gaussian')
        residual=(mag-profile(best,ph,'gaussian'))/np.hypot(err,.012)
        null=best.copy()
        midpoint=.5*(best[2]+best[5]);null[2],null[5]=midpoint-.25,midpoint+.25
        baseline=profile(null,ph,'gaussian')
        trials=[]
        for _ in range(a.draws):
            noise=rng.choice(residual,size=len(err))*np.hypot(err,.012)
            trials.append(derived(fit(ph,baseline+noise,err,null,'gaussian'))[2])
        x.append((best[5]-best[2],np.array(trials)))
    observed=x[1][0]-x[0][0];delta=x[1][1]-x[0][1]
    a.output.write_text(json.dumps({'draws':a.draws,'seed':803122,'null_separation':.5,
        'description':'Actual OGLE-IV early and late timestamps, each era retains its own fitted widths/depths/baseline and mean eclipse-center phase, with primary-to-secondary spacing fixed at 0.5 cycles. Add independent draws from that era’s empirical standardized residuals. Refit both centers freely. This single-target noise-and-cadence check does not calibrate a survey false alarm rate or mimic phase-dependent residual correlations.',
        'observed_delta':observed,'central95_null_delta':np.quantile(delta,[.025,.975]).tolist(),
        'largest_abs_null_delta':float(np.max(abs(delta))),
        'n_at_least_observed_absolute':int((abs(delta)>=abs(observed)).sum())},indent=2)+'\n')

if __name__=='__main__':main()
