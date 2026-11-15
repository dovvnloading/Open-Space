"""Adversarial check for a sampled high-order harmonic of the primary rotation."""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from measure import residuals, event_rows


def residual_at_frequency(t, y, f, harmonics):
    phase = 2*np.pi*t[:,None]*f*np.arange(1,harmonics+1)
    x = t / np.ptp(t) - .5
    design = np.column_stack([np.ones(len(t)),x,x*x,np.sin(phase),np.cos(phase)])
    return y-design @ np.linalg.lstsq(design,y,rcond=None)[0]


def main(release, out, destination):
    meta = {r[0]:float(r[3]) for line in (release/'release.merge').read_text().splitlines() if (r:=line.split())}
    with (out/'measurements.csv').open() as f:
        m = {int(r['mpc']):r for r in csv.DictReader(f)}
    rows=[]
    for i in [29489,3177]:
        r=m[i]; a,_,_=residuals(release,i,meta[str(i)])
        jd=a[:,0]; t=jd-jd[0]; span=np.ptp(t)
        mag=a[:,3]-5*np.log10(a[:,8]*a[:,9])
        f0=24/float(r['spin_h']); p=float(r['recurrence_h']); center=float(r['first_center_jd'])
        w=float(r['box_width_h'])/p
        # Native prime-mission sampling is 48 per day; archive JD is rounded.
        orders=np.arange(1,33)
        aliases=np.abs(orders*f0-48*np.rint(orders*f0/48))
        closest=np.argmin(abs(aliases-24/p))
        for harmonics in [6,14,28]:
            # Scan first to choose a local basin, then refine. Same narrow
            # search bounds for each order; high orders are deliberately flexible.
            freqs=np.linspace(f0-.0025,f0+.0025,51)
            def loss(f):
                err=residual_at_frequency(t,mag,f,harmonics)
                return np.mean(err**2)
            scores=np.array([loss(freq) for freq in freqs])
            j=int(np.argmin(scores))
            bounds=(freqs[max(0,j-1)],freqs[min(j+1,len(freqs)-1)])
            opt=minimize_scalar(loss,bounds=bounds,method='bounded',options={'xatol':1e-9})
            f=opt.x if opt.fun <= scores[j] else freqs[j]
            resid=residual_at_frequency(t,mag,f,harmonics)
            ev=event_rows(jd-center,resid,p,0,w)
            rows.append(dict(mpc=i,spin_harmonics=harmonics,fitted_spin_cpd=f,
                recurrence_cpd=24/p,rayleigh_resolution_cpd=1/span,
                nearest_harmonic_le32=int(orders[closest]),near_alias_cpd=aliases[closest],
                alias_mismatch_cpd=abs(aliases[closest]-24/p),
                fit_rms_mag=np.std(resid),median_depth_mag=np.median(ev[:,3]),
                n_positive=int((ev[:,3]>0).sum()),n_events=len(ev)))
    destination.mkdir(parents=True,exist_ok=True)
    with (destination/'spin_alias.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader()
        writer.writerows([{k:round(v,8) if isinstance(v,(float,np.floating)) else v for k,v in row.items()} for row in rows])


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('release',type=Path)
    p.add_argument('results',type=Path)
    p.add_argument("--out",type=Path,default=Path(__file__).resolve().parent)
    args=p.parse_args()
    main(args.release,args.results,args.out)
