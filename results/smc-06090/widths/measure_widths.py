"""Test whether the existing TESS aperture supports a duration asymmetry.

The metric is the fitted difference in full width at half maximum in flux,
in minutes (secondary minus primary). Fit independent centers and depths.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'joint'))
from tess_template import fit, GROUPS
PERIOD=1.1696719
TO_MINUTES=2*np.sqrt(2*np.log(2))*PERIOD*1440

def blocks(cycles):
    runs=np.split(np.arange(len(cycles)),np.flatnonzero(np.diff(cycles)!=1)+1)
    return [run[j:j+3] for run in runs for j in range(0,len(run),3)]

def measurement(S,N):
    fitted=fit(S,N,wide=True)
    a,b=fitted[-2:]*TO_MINUTES
    return np.array([a,b,b-a])

def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--draws',type=int,default=600)
    args=cli.parse_args()
    path=HERE.parent/'joint'/'inputs'
    meta=json.loads((path/'tess_inputs.json').read_text())
    rng=np.random.default_rng(4167904)
    out={}
    with np.load(path/'orbital_folds.npz') as data:
        for group,sectors in GROUPS.items():
            S,N,groups=[],[],[]
            for s in sectors:
                n=data[f's{s}_counts'];dt=meta[str(s)]['cadence_days']
                ok=(n.sum(axis=1)*dt>.65*PERIOD)&((n>0).sum(axis=1)>=min(55,int(.65*PERIOD/dt)))
                S.append(data[f's{s}_sums'][ok]);N.append(n[ok])
                groups.append(blocks(data[f's{s}_cycles'][ok]))
            aperture={size:measurement([a[:,:,i] for a in S],N).tolist()
                      for i,size in enumerate(['3x3','four','one'])}
            primary=[a[:,:,0] for a in S]
            draws=[]
            for i in range(args.draws):
                ix=[np.concatenate([g[j] for j in rng.integers(0,len(g),size=len(g))]) for g in groups]
                draws.append(measurement([s[k] for s,k in zip(primary,ix)],[n[k] for n,k in zip(N,ix)]))
            bounds=np.percentile(draws,[2.5,97.5],axis=0)
            out[group]=dict(apertures_minutes=aperture,ci95_minutes=bounds.T.tolist(),n_bootstrap=args.draws)
            print(group,aperture['3x3'],bounds[:,2],flush=True)
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':
    main()
