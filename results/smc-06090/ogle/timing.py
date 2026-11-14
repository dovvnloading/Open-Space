"""Reproduce the eclipse timing reanalysis using primary OGLE photometry.

Command: OPENBLAS_NUM_THREADS=1 python timing.py --output measurements.json
Dependencies: NumPy and SciPy. Source input files are downloaded into the cache.
"""
import argparse
import csv
import gzip
import hashlib
import json
import shutil
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.optimize import least_squares

OBJECT = "OGLE-SMC-ECL-06090"
BASE = "https://www.astrouw.edu.pl/ogle/"
REPO_REV = "5f8bc79e188ff1940d87a59e43fe2e96a9083059"
FILES = {
    "ogle3.dat": (BASE+"ogle3/OIII-CVS/smc/ecl/phot/I/OGLE-SMC-ECL-6090.dat", "d57614117ef140adf1413baab42f92182b9c39752eed39706b6fc040e203a897"),
    "ogle4.dat": (BASE+"ogle4/OCVS/smc/ecl/phot/I/"+OBJECT+".dat", "429caaf74f2086f3a8f53c71a91bfb3c07049dbc2e81122c06646387c8fa21a0"),
    "ecl.dat": (BASE+"ogle4/OCVS/smc/ecl/ecl.dat", "df02fb6ac66e3377e326f1c1b54e7b20230e15952540c3e7ad16f46f42073185"),
    "screen.csv.gz": ("https://raw.githubusercontent.com/dovvnloading/Open-Space/"+REPO_REV+"/results/smc-eclipse-depths/screen.csv.gz", "e4a384db50b1be35d4a35ef54baadaff3370e80e8ea8cb71b7de4fce0a8466e6")
}
SPLITS = [
    ("OGLE-III early", "ogle3.dat", 0, 2453500),
    ("OGLE-III late", "ogle3.dat", 2453500, np.inf),
    ("OGLE-IV 2010–14", "ogle4.dat", 0, 2457000),
    ("OGLE-IV 2014–18", "ogle4.dat", 2457000, 2458300),
    ("OGLE-IV 2018–22", "ogle4.dat", 2458300, 2459700),
    ("OGLE-IV 2022–24", "ogle4.dat", 2459700, np.inf)
]
def obtain(cache):
    cache.mkdir(parents=True, exist_ok=True)
    for name, (url, expected) in FILES.items():
        path = cache / name
        if not path.exists():
            root = Path(__file__).resolve().parent
            local = (root.parent / 'joint' / 'inputs' / name if name != 'screen.csv.gz'
                     else root.parent.parent / 'smc-eclipse-depths' / name)
            if local.is_file():
                shutil.copyfile(local, path)
            else:
                with urlopen(url, timeout=100) as response:
                    path.write_bytes(response.read())
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"input hash mismatch for {path}")

def select(cache):
    with gzip.open(cache / 'screen.csv.gz','rt') as f: rows = list(csv.DictReader(f))
    opposite=[]; cuts=[]
    for row in rows:
        a,b = [float(row[col]) for col in ('primary_delta_mag','secondary_delta_mag')]
        if a*b>=0: continue
        opposite.append(row)
        if min(abs(a),abs(b))>=.1:
            score = min(abs(a)/float(row['primary_nominal_error']),abs(b)/float(row['secondary_nominal_error']))
            cuts.append((score,row['id']))
    cuts=sorted(cuts, reverse=True)
    assert cuts[0][1]==OBJECT
    return {'screen_n':len(rows), 'opposite_sign_n':len(opposite), 'large_opposite_ranked':cuts}

def signed_phase(x): return (x+.5)%1-.5
def profile(par, ph, shape):
    level,d1,c1,w1,d2,c2,w2=par
    def dip(c,w):
        u=signed_phase(ph-c)/w
        if shape=='gaussian': return np.exp(-.5*u**2)
        return np.maximum(1-abs(u)/2.5,0)
    return level+d1*dip(c1,w1)+d2*dip(c2,w2)

def fit(ph,mag,err,start,shape):
    baseline=np.quantile(mag,.25)
    return least_squares(lambda par:(profile(par,ph,shape)-mag)/np.hypot(err,.012),
           start, bounds=([baseline-.25,.005,-.12,.008,.005,.38,.008],
                         [baseline+.25,1.5,.12,.09,1.5,.62,.09]),
           loss='soft_l1',f_scale=1.8, max_nfev=120).x

def derived(par):
    return np.array([par[2],par[5],par[5]-par[2],par[1],par[4],par[0]])

def summary(t, ph, mag, err, rng, draws):
    nights, labels=np.unique(np.floor(t), return_inverse=True)
    groups=[np.flatnonzero(labels==k) for k in range(len(nights))]
    initial=[np.quantile(mag,.25),.4,.0,.023,.4,.51,.023]
    par=fit(ph,mag,err,initial,'gaussian')
    tri=fit(ph,mag,err,par,'triangular')
    samples=[]
    for i in range(draws):
        ids=np.concatenate([groups[k] for k in rng.integers(0,len(nights),len(nights))])
        samples.append(derived(fit(ph[ids],mag[ids],err[ids],par,'gaussian')))
    dist=np.asarray(samples)
    return {'used_rows':len(t),'nights':len(nights),'hjd_range':[t.min(),t.max()],
            'gaussian':derived(par).tolist(),'night_bootstrap_95':np.quantile(dist,[.025,.975],axis=0).T.tolist(),
            'triangle':derived(tri).tolist()}, dist

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache',default='timing-inputs', type=Path)
    ap.add_argument('--draws',default=800, type=int)
    ap.add_argument('--seed',default=614823,type=int)
    ap.add_argument('--output',default='measurements.json',type=Path)
    args=ap.parse_args();obtain(args.cache)
    row=next(l.split() for l in (args.cache / 'ecl.dat').read_text().splitlines() if l.split()[0]==OBJECT)
    period,epoch=map(float,row[3:5])
    measurements={'object':OBJECT, 'catalog_period_d':period,'catalog_epoch_hjd':epoch,
        'selection':select(args.cache),'draws':args.draws,'seed':args.seed,
        'columns':['primary_phase','secondary_phase','separation_cycles','primary_depth_mag','secondary_depth_mag','baseline_mag'],
        'inputs':FILES, 'epochs':{},'input_rows_usable_total':{},'comparisons':{}}
    rng=np.random.default_rng(args.seed);sims=[];centers=[]
    for name,source,low,high in SPLITS:
        t,m,e=np.loadtxt(args.cache/source,usecols=(0,1,2)).T
        if source=='ogle3.dat':t+=2450000
        ok = np.isfinite(t+m+e)&(e>0)&(e<=.1)
        measurements['input_rows_usable_total'][source]=[int(ok.sum()), len(t)]
        keep=ok&(t>=low)&(t<high)
        times=t[keep]; ph=signed_phase((times-epoch)/period)
        rec,sample=summary(times,ph,m[keep],e[keep],rng,args.draws)
        sims.append(sample);centers.append(np.asarray(rec['gaussian']))
        measurements['epochs'][name]=rec
    for label,i,j in [('OGLE-IV last minus first',2,5),('OGLE-III last minus first',0,1),('OGLE-IV first minus OGLE-III first',0,2)]:
        diff=sims[j]-sims[i]
        measurements['comparisons'][label]={'delta':(centers[j]-centers[i]).tolist(),
            'night_bootstrap_95':np.quantile(diff,[.025,.975],axis=0).T.tolist()}
    args.output.write_text(json.dumps(measurements,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':
    main()
