"""Measure four folded-light-curve intervals for the entire catalog-defined population."""
import argparse
import csv
import hashlib
import json
import tarfile
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
TARGET = 'OGLE-SMC-ECL-06090'
BASE = {3: 'https://www.astrouw.edu.pl/ogle/ogle3/OIII-CVS/smc/ecl/',
        4: 'https://www.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/'}
WINDOWS = [(3, 0, 2453500), (3, 2453500, 2455200),
           (4, 0, 2457000), (4, 2459700, 2500000)]


def catalog(cache):
    old = {}
    for line in (cache / 'ogle3-ident.dat').read_text().splitlines():
        name, field, local = line.split()[:3]
        old[f'{field}.{local}'] = name
    identity = {s[0]: s for s in map(str.split, (cache / 'ogle4-ident.dat').read_text().splitlines())}
    stars = []
    for row in (cache / 'ogle4-ecl.dat').read_text().splitlines():
        name, i, _, period, epoch, d1, d2 = row.split()
        try:
            mag, p, t0, depth = float(i), float(period), float(epoch), min(float(d1),float(d2))
        except ValueError:
            continue
        ident = identity[name]
        if ident[1]=='NC' and abs(mag-18.225)<=.75 and .75*1.1696719<=p<=1.5*1.1696719 and depth>=.15:
            prior = next((old[v] for v in ident[5:] if v in old), None)
            stars.append((name, prior, p, t0))
    return stars


def materialize(cache, stars):
    for s in (4,3):
        names = {row[0 if s==4 else 1] for row in stars}
        target_dir = cache / str(s)
        target_dir.mkdir(exist_ok=True)
        with tarfile.open(cache/f'ogle{s}-phot.tar.gz','r|gz') as archive:
            for member in archive:
                path = Path(member.name)
                if member.isfile() and path.parent.as_posix()=='phot/I' and path.stem in names:
                    (target_dir/path.name).write_bytes(archive.extractfile(member).read())


def prepare(cache):
    cache.mkdir(exist_ok=True,parents=True)
    for era, fname in ((3,'ident.dat'),(4,'ident.dat'),(4,'ecl.dat'),(3,'phot.tar.gz'),(4,'phot.tar.gz')):
        path=cache/f'ogle{era}-{fname}'
        if not path.exists():
            urlretrieve(BASE[era]+fname,path)
    stars=catalog(cache)
    materialize(cache,stars)
    receipt = {p.name:hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
               for p in sorted(cache.glob('ogle*')) if p.is_file()}
    (HERE/'sources.json').write_text(json.dumps({'surveys':BASE,'sha256':receipt},indent=2)+'\n')


def wrap(x):
    return (x+.5)%1-.5


def measure(data, period, epoch):
    t, m, errors = data.T
    n=len(m)
    if n<40:
        return {'n':n,'quality':'few points'}
    phase=wrap((t-epoch)/period)
    harmonic=np.cos(4*np.pi*phase)
    sigma=np.hypot(errors,.012)
    def evaluate(x,jacobian=False):
        base,ell,d1,c1,w1,d2,c2,w2=x
        u1,u2=wrap(phase-c1)/w1,wrap(phase-c2)/w2
        g1,g2=np.exp(-.5*u1**2),np.exp(-.5*u2**2)
        if jacobian:
            return np.column_stack((np.ones(n),harmonic,g1,d1*g1*u1/w1,d1*g1*u1**2/w1,
                                    g2,d2*g2*u2/w2,d2*g2*u2**2/w2))/sigma[:,None]
        return (base+ell*harmonic+d1*g1+d2*g2-m)/sigma
    base=np.percentile(m,35)
    lo=np.array([base-.4,-.20,.001,-.21,.008,.001,.29,.008])
    hi=np.array([base+.4,.20,1.8,.21,.12,1.8,.71,.12])
    initial=[(0,.5)]
    # Independent coarse template locations reduce local minima for drifting catalog ephemerides.
    def locate(center):
        bins=center+np.arange(-.18,.181,.015)
        dist=wrap(phase[:,None]-bins)
        w=np.exp(-.5*(dist/.025)**2)/sigma[:,None]**2
        depths=(w*(np.minimum(m,base+1.5)-base)[:,None]).sum(axis=0)/np.maximum(w.sum(axis=0),1e-10)
        return bins[np.argmax(depths)]
    alt=(locate(0.),locate(.5))
    if max(abs(alt[0]),abs(alt[1]-.5))>.022:
        initial.append(alt)
    best=None
    for c1,c2 in initial:
        x=[base,.015,.28,c1,.03,.23,c2,.03]
        res=least_squares(evaluate,x,jac=lambda a:evaluate(a,True),bounds=(lo,hi),
                          loss='soft_l1',f_scale=1.8,max_nfev=180)
        if best is None or res.cost<best.cost:
            best=res
    pars=best.x
    separation=float(pars[6]-pars[3])
    nights=[len(set(np.floor(t[abs(wrap(phase-c))<.06]))) for c in pars[[3,6]]]
    fit_good=(best.success and min(pars[[2,5]])>=.05 and
              np.all(np.minimum(pars-lo,hi-pars)[[2,3,4,5,6,7]]>1.e-4))
    jac=best.jac
    cov=np.linalg.pinv(jac.T@jac)
    v=np.zeros(8);v[3]=-1;v[6]=1
    scale=max(1,1.4826*np.median(abs(best.fun-np.median(best.fun))))
    se=float(np.sqrt(max(0,v@cov@v))*scale)
    quality='usable' if fit_good and np.linalg.matrix_rank(jac)==8 and min(nights)>=5 else 'shape/coverage'
    return dict(n=n,quality=quality,sep=separation,se=se,core1=nights[0],core2=nights[1],cost=best.cost)


def one_star(star, cache):
    name, prior, p, epoch=star
    rows=[]
    for survey, name_in_survey in ((3,prior),(4,name)):
        path=cache/str(survey)/f'{name_in_survey}.dat'
        if name_in_survey and path.exists():
            data=np.loadtxt(path,usecols=(0,1,2),ndmin=2)
            data=data[np.isfinite(data).all(axis=1)&(data[:,2]>0)&(data[:,2]<=.10)]
            if survey==3: data[:,0]+=2450000
        else:
            data=np.empty((0,3))
        for s,a,b in WINDOWS:
            if s==survey:
                values=data[(data[:,0]>=a)&(data[:,0]<b)]
                rows.append(dict(name=name,prior=prior or '',period=p,survey=s,from_hjd=a,to_hjd=b,
                                 **measure(values,p,epoch)))
    return rows


def summarize(rows,limit):
    stars={}
    for row in rows:
        stars.setdefault(row['name'],[]).append(row)
    def diff(ends):
        if all(r['quality']=='usable' and float(r['se'])<=limit for r in ends):
            return float(ends[1]['sep'])-float(ends[0]['sep'])
    usable={k:(diff(x[:2]),diff(x[2:])) for k,x in stars.items()}
    first,last=usable[TARGET]
    good=[(k,y) for k,(x,y) in usable.items() if y is not None]
    both=[(k,x,y) for k,(x,y) in usable.items() if x is not None and y is not None]
    return dict(precision_cut=limit, target=[first,last], retained_IV=len(good),retained_both=len(both),
                iv_at_least=[k for k,y in good if abs(y)>=abs(last)-1e-12] if last is not None else None,
                reversal_at_least=[k for k,x,y in both if x*y<0 and abs(x)>=abs(first)-1e-12 and abs(y)>=abs(last)-1e-12] if first is not None and last is not None else None,
                comparison_IV_abs_quantiles=np.percentile([abs(y) for k,y in good if k!=TARGET],[50,90,95,99]).tolist())


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--cache',type=Path,default=HERE.parent/'cache')
    cli.add_argument('--prepare',action='store_true')
    cli.add_argument('--workers',type=int,default=3)
    args=cli.parse_args()
    if args.prepare:
        prepare(args.cache)
    with ProcessPoolExecutor(args.workers) as pool:
        rows=[r for block in pool.map(partial(one_star,cache=args.cache),catalog(args.cache),chunksize=4) for r in block]
    cols='name prior period survey from_hjd to_hjd n quality sep se core1 core2 cost'.split()
    with (HERE/'epoch_measurements.csv').open('w') as out:
        writer=csv.DictWriter(out,fieldnames=cols);writer.writeheader();writer.writerows(rows)
    result={str(limit):summarize(rows,limit) for limit in [.0075,.01,.015,.020,.030]}
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['0.015'],flush=True)


if __name__=='__main__':
    main()
