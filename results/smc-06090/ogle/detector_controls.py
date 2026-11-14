"""Membership and timing checks for the target's OGLE-IV chip.

Run from any directory: python detector_controls.py --cache PATH --out output.json.
The input files are publicly served by the OGLE collaboration.
"""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.optimize import least_squares

BASE='https://www.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/'
TARGET='OGLE-SMC-ECL-06090'
INTERVALS=((0,2457000),(2459700,1e9))


def get(name,cache):
    path=cache/name
    if not path.exists():
        path.parent.mkdir(parents=True,exist_ok=True)
        with urlopen(BASE+name,timeout=60) as response:
            path.write_bytes(response.read())
    return path


def objparams(catalog,ident):
    rows={}
    tags={s[:18]:(s[20:23].strip(),s[49:65].strip().rsplit('.',1)[0]) for s in ident.splitlines()}
    for row in catalog.splitlines():
        name,i,_,p,t,d1,d2=row.split()
        if '-' not in (i,p,t,d1,d2):
            rows[name]=[float(v) for v in (i,p,t,d1,d2)]
    it,pt,_,_,_=rows[TARGET]
    return {n:row for n,row in rows.items() if
        tags[n]==tags[TARGET] and abs(row[0]-it)<=1 and 0.5*pt<=row[1]<=2*pt
        and min(row[-2:])>=0.10}


def measurement(data,lo,hi,p,t0):
    t,m,e=data.T
    good=(t>=lo)&(t<hi)&np.isfinite(t+m+e)&(e>0)&(e<=.1)
    t,m,e=t[good],m[good],e[good]
    ph=((t-t0)/p+.25)%1-.25
    if len(t)<25:
        return {'n':len(t),'valid':False}
    def mdl(par):
        b,c=par[:2]; y=b+c*np.cos(4*np.pi*ph)
        for a,mu,s in (par[2:5],par[5:8]):
            u=((ph-mu+.5)%1-.5)/s
            y=y+a*np.exp(-.5*u**2)
        return y
    b=np.median(m)
    fit=least_squares(lambda pars:(mdl(pars)-m)/np.maximum(e,.012),
                      [b,0,.3,0,.03,.25,.5,.03],
                      bounds=([b-.3,-.15,0,-.18,.009,0,.32,.009],
                              [b+.3,.15,1.5,.18,.12,1.5,.68,.12]),loss='soft_l1')
    par=fit.x
    nights=[len(np.unique(np.floor(t[np.abs((ph-mu+.5)%1-.5)<.06]))) for mu in par[[3,6]]]
    return {'n':len(t),'core_nights':nights, 'fits':par.tolist(),'separation':par[6]-par[3],
            'valid':bool(min(nights)>=5 and min(par[[2,5]])>=.05 and
                         min(par[[4,7]])>.00901 and max(par[[4,7]])<.11999 and
                         abs(par[3])<.1799 and abs(par[6]-.5)<.1799)}


def main(cache, output):
    cat=get('ecl.dat',cache); ident=get('ident.dat',cache)
    hashes={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in (cat,ident)}
    result={}
    for name,(i,p,t,d1,d2) in objparams(cat.read_text(),ident.read_text()).items():
        path=get('phot/I/'+name+'.dat',cache)
        hashes['phot/I/'+path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        data=np.loadtxt(path)
        eras=[measurement(data,*lim,p,t) for lim in INTERVALS]
        diff=eras[1]['separation']-eras[0]['separation'] if all(a['valid'] for a in eras) else None
        result[name]={'I':i,'period':p,'early_late':eras,'delta_separation':diff}
    output.write_text(json.dumps({'base':BASE,'hashes':hashes,'objects':result},indent=2)+'\n')
    print([(k,r['delta_separation']) for k,r in result.items()])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,default=Path('controls_cache'))
    p.add_argument('--out',type=Path,default=Path('detector_controls.json'))
    a=p.parse_args()
    main(a.cache,a.out)
