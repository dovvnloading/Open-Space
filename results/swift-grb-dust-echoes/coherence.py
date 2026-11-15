"""Rank the whole screen by the weakest of two time and four angular partitions."""
import argparse
from pathlib import Path

import numpy as np
from scipy.special import xlogy
from search import rows, write_csv

GRID = np.exp(np.arange(np.log(30),np.log(16000),.025))


def scan(a):
    t,x,y,e = a
    keep = (e>=70)&(e<=300)
    if not keep.any():
        return (0.,0.,0)
    epoch = np.argsort(t)
    snap = np.empty(len(t),dtype=int)
    snap[epoch] = np.r_[0,np.cumsum(np.diff(t[epoch])>500)]
    t,x,y,snap = t[keep],x[keep],y[keep],snap[keep]
    radius2=x*x+y*y
    angle=np.floor((np.arctan2(y,x)+np.pi)*2/np.pi).astype(int)%4
    # At every distance, search the completely contained rings in six partitions.
    # Min of counts-only log-likelihood scores disfavors a single stationary source.
    best = (0., 0., 0)
    for dist in GRID:
        radii2=826.6*t/dist
        covered=(.7*radii2>2500)&(1.3*radii2<550**2)
        ratio=radius2/radii2
        allband=covered&(ratio>=.7)&(ratio<1.3)
        on=allband&(ratio>=.9)&(ratio<1.1)
        n,m=(np.r_[np.bincount(snap[k]%2,minlength=2),
                   np.bincount(angle[k],minlength=4)] for k in (on,allband))
        if np.any(m==0) or np.any(n<=m/3):
            continue
        scores=2*(xlogy(n,3*n/m)+xlogy(m-n,1.5*(m-n)/m))
        best=max(best,(float(np.sqrt(np.min(scores))),float(dist),int(on.sum())))
    return best


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("photons",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    output=[]
    for row in rows("sample.csv"):
        path=args.photons/(row["name"]+".npy")
        value=scan(np.load(path))
        output.append(dict(name=row["name"],min_partition_score=round(value[0],4),
                           distance_pc=round(value[1],2),on=value[2]))
    write_csv(args.output,sorted(output,key=lambda r:-r["min_partition_score"]))


if __name__=="__main__":
    main()
