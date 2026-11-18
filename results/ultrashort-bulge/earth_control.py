"""Recover the published OGLE-2016-BLG-1928 signal outside the search years."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.optimize import least_squares
import fit
from replay_guard import fresh,external

ROOT=Path(__file__).resolve().parent
URL='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J%2FApJ%2F903%2FL11&-out.max=unlimited&-out.all=1'


def run(raw,library_path=None):
    data=[]
    for line in raw.decode().splitlines():
        cols=line.strip().split('\t')
        if len(cols)>=5 and cols[0].strip().isdigit():
            data.append((float(cols[1]),float(cols[2]),float(cols[3]),cols[4].strip()))
    t,m,e,site=map(np.asarray,zip(*data))
    f,se=fit.lens.to_flux(m,e)
    g=np.unique(site,return_inverse=True)[1]
    fit.setup(library_path)
    # Use the published observing night as a position hint and no timescale hint.
    summary=fit.screen(t,f,se,g,2457557.8,np.nan,2016,False)
    vb=fit.engine()
    # Hold-out validation for the isolated/no-blend source hypothesis; baseline
    # determined from exposures away from the published event. Distinct streams
    # have independently measured baselines. Free blending is tested separately.
    ref=2457557.8;days=t-ref; near=abs(days)<1
    baseline=np.array([np.median(f[(g==j)&~near]) for j in np.unique(g)])
    tn=days[near];fn=f[near];en=se[near];gn=g[near]
    def fit_one(finite,blend):
        def err(x):
            dt,lu,lte=x[:3];u,te=np.exp([lu,lte])
            if finite:amp=np.array(vb.ESPLLightCurve([lu,lte,dt,x[3]],tn.tolist())[0])
            else:amp=fit.lens.magnification(u*u+((tn-dt)/te)**2)
            fraction=x[-len(baseline):] if blend else np.ones(len(baseline))
            return (baseline[gn]*(1+fraction[gn]*(amp-1))-fn)/en
        opts=[]
        for te0 in [.025,.07,.18]:
            for u0 in [.15,.9]:
                x=[0.,np.log(u0),np.log(te0)]+([np.log(3.)] if finite else [])
                lo=[-.25,-6,np.log(.005)]+([np.log(.02)] if finite else [])
                hi=[.25,np.log(9),np.log(5)]+([np.log(20)] if finite else [])
                if blend:
                    x+=list(np.ones(len(baseline)));lo+=[.01]*len(baseline);hi+=[1.15]*len(baseline)
                opts.append(least_squares(err,x,bounds=(lo,hi),max_nfev=400,diff_step=1e-5))
        opt=min(opts,key=lambda s:s.cost);v=opt.x
        return dict(t0=ref+v[0],u0=np.exp(v[1]),te_d=np.exp(v[2]),rho=np.exp(v[3]) if finite else 0.,
                    source_fraction=v[-len(baseline):].tolist() if blend else [1.]*len(baseline),
                    chi2=2*opt.cost)
    return {'event':'OGLE-2016-BLG-1928','doi':'10.3847/2041-8213/abbfad','URL':URL,
            'response_sha256':hashlib.sha256(raw).hexdigest(),'input_rows':len(data),
            'sites':np.unique(site).tolist(),'observations_within_1d':int(near.sum()),
            'grid':summary,'point_no_blend':fit_one(False,False),'finite_no_blend':fit_one(True,False),
            'point_free_blend':fit_one(False,True),'finite_free_blend':fit_one(True,True),
            'published':{'te_d':.0288,'rho':3.39,'theta_e_uas':.842,'theta_e_err_uas':.064}}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--grid-library',type=Path,required=True)
    args=parser.parse_args()
    input_path=external(args.input);raw=input_path.read_bytes()
    if raw[:2]==b'\x1f\x8b':raw=gzip.decompress(raw)
    expected='819169bc4c6ab82941c2c01b48c8ba9c331bd55735348d30f2d5715bc2cb8c0c'
    if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('frozen positive-control input mismatch')
    out=fresh(args.out)
    result=run(raw,external(args.grid_library))
    (out/'earth-recovery.json').write_text(json.dumps(result,indent=2,default=lambda v:v.item())+'\n')
