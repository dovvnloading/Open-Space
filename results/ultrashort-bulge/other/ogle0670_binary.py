"""Transfer converged OGLE-0670 binary geometry across flux ceilings and bin choices.

Hold lens shape fixed at the common 30-start solution and independently solve for
two flux terms. Compare with continuously reoptimized FSPL under each ceiling.
"""
import csv
import json

import numpy as np
from ogle0670 import NAME,OUT,read_curve,write_csv
from pathlib import Path
import fit

if __name__=='__main__':
    shape=json.loads((OUT/'ogle0670-models.json').read_text())
    p=Path('.ogle-cache/raw/2026/blg-0670.dat')
    t,y,e,g,diff=read_curve((NAME,p,2461193.333,8.497))
    t,y,e,g=fit.clean(t,y,e,g);t-=shape['reference_jd']
    q=shape['models']['bin']['vector']
    v=fit.engine();v.Tol=.001
    near=abs(t)<=30;far=abs(t)>90
    z=(y[far]-np.median(y[far]))/e[far];escale=max(1.,np.median(abs(z))/.67448975)
    e*=escale;ok=abs(z)<6*escale;w0=1/e[far][ok]**2
    f0=np.average(y[far][ok],weights=w0)
    rows=[]
    for r in csv.DictReader((OUT/'ogle0670_profiles.csv').open()):
        d=float(r['bin_d']);factor=(np.inf if r['baseline_flux_factor']=='unconstrained' else float(r['baseline_flux_factor']))
        tn,yn,en=t[near],y[near],e[near]
        if d:
            g=np.unique(np.floor(tn/d),return_inverse=True)[1]
            weights=1/en**2;tot=np.bincount(g,weights=weights)
            tn=np.bincount(g,weights=tn*weights)/tot;yn=np.bincount(g,weights=yn*weights)/tot;en=1/np.sqrt(tot)
        a=np.asarray(v.BinaryLightCurve([q[4],q[5],np.exp(q[2]),q[6],q[3],q[1],q[0]],tn.tolist())[0])-1
        a=np.r_[a,0.];wc=np.r_[1/en**2,np.sum(w0)];fy=np.r_[yn,f0];fy-=np.average(fy,weights=wc);a-=np.average(a,weights=wc)
        fs=np.clip(np.dot(wc*a,fy)/np.dot(wc,a*a),0.,factor*np.median(y))
        chi2=float(np.sum(wc*(fy-fs*a)**2))
        rows.append(dict(event=NAME,baseline_flux_factor=r['baseline_flux_factor'],bin_d=d,
                         fixed_binary_chi2=round(chi2,3),binary_source_fraction=round(fs/np.median(y),6),
                         profiled_fs_chi2=r['fs_chi2'],fs_minus_fixed_bin=round(float(r['fs_chi2'])-chi2,3)))
    write_csv(OUT/'ogle0670_fixed_binary.csv',rows)
