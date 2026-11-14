"""Screen bright short-period LMC binaries, reserving 2022-2024 as a time holdout.

Run with python screen.py --cache PATH. See README.md for source verification.
"""
from pathlib import Path
import argparse
import csv
import tarfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT=Path(__file__).resolve().parent
RAW=ROOT/'raw'
DATA=ROOT/'selected'
EARLY=(2455197.5,2457023.5)  # 2010–2014
TRAIN=(2457754.5,2458850.0)  # 2017–2019
HOLD=(2459580.5,2460676.5)   # 2022–2024


def catalog():
    ident={s[:18]:s[20:23].strip() for s in (RAW/'ident.dat').read_text().splitlines()}
    chosen={}
    for line in (RAW/'ecl.dat').read_text().splitlines():
        name,i,v,per,epoch,d1,d2=line.split()
        if '-' in (i,v,per,epoch,d1,d2):
            continue
        if (ident[name]=='NC' and float(i)<=18.5 and 0.7<float(per)<6
                and float(d1)>=0.16 and float(d2)>=0.12):
            chosen[name]=(float(per),float(epoch))
    return chosen


def extract(chosen):
    DATA.mkdir(exist_ok=True)
    if (DATA/'complete').exists():
        return
    with tarfile.open(RAW/'phot.tar.gz','r|gz') as tar:
        for member in tar:
            path=Path(member.name)
            if member.isfile() and path.parent.as_posix()=='phot/I' and path.stem in chosen:
                (DATA/path.name).write_bytes(tar.extractfile(member).read())
    (DATA/'complete').touch()


def per_era(data, limits, p, ep):
    t,m,e=data.T
    ok=(t>=limits[0])&(t<limits[1])&np.isfinite(m)&(e>0)&(e<=0.10)
    t,m,e=t[ok],m[ok],e[ok]
    ph=((t-ep)/p+0.5)%1.-0.5
    q=(abs(ph)>.16)&(abs(ph)<.34)
    if q.sum()<20:
        return
    # Equal weight per night in each window. Floor HJD gives a consistent
    # midday boundary for southern-night OGLE measurements.
    def medians(idx):
        n=np.floor(t[idx]).astype(int)
        if len(np.unique(n))==len(n):
            return m[idx]
        return np.array([np.median(m[idx][n==night]) for night in np.unique(n)])
    stats=[]
    for window in (q,abs(ph)<.03,abs(ph)>.47):
        vals=medians(window)
        if len(vals)<8:
            return
        value=np.median(vals)
        err=1.2533*max(1.4826*np.median(abs(vals-value)),np.median(e[window]),.008)/np.sqrt(len(vals))
        stats.append((value,err,len(vals)))
    return stats


def screen_pair(item):
    name,(period,epoch)=item
    path=DATA/(name+'.dat')
    if not path.exists():
        return
    data=np.loadtxt(path,ndmin=2)
    eras=[per_era(data,lim,period,epoch) for lim in (EARLY,TRAIN)]
    if None in eras:
        return
    a,b=eras
    delta=[(b[j][0]-b[0][0])-(a[j][0]-a[0][0]) for j in (1,2)]
    z=[delta[j-1]/np.sqrt(sum(x[k][1]**2 for x in (a,b) for k in (0,j))) for j in (1,2)]
    score=min(map(abs,z)) if delta[0]*delta[1]>0 else 0.
    return [name,period,epoch,score,*delta,*z,b[0][0]-a[0][0]]


if __name__=='__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--cache', type=Path, required=True,
                     help='Directory with raw/ and writable selected/ and screen.csv')
    args = cli.parse_args()
    RAW = args.cache / 'raw'
    DATA = args.cache / 'selected'
    chosen=catalog()
    print('Eligible',len(chosen),flush=True)
    extract(chosen)
    with ProcessPoolExecutor(max_workers=5) as ex:
        rows=[r for r in ex.map(screen_pair,chosen.items(),chunksize=20) if r]
    rows.sort(key=lambda r:(-r[3],r[0]))
    with (args.cache/'screen.csv').open('w') as f:
        w=csv.writer(f)
        w.writerow(['star','period','epoch','rank_score','change_primary','change_secondary',
                    'z_primary','z_secondary','change_quadrature'])
        w.writerows(rows)
    print('With files',len(list(DATA.glob('*.dat'))),'adequate',len(rows))
    remarks=dict(line.split(' ',1) for line in (RAW/'remarks.txt').read_text().splitlines() if ' ' in line)
    for r in rows[:25]: print(r,remarks.get(r[0],''))
