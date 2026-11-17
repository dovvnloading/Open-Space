"""Exact duration with two distinct exposures on each side of a uniform trial peak."""
import argparse
import csv
import json
from pathlib import Path
import multiprocessing as mp
import numpy as np
from run_grid import inputs, read_curve
from fit import clean

WIDTHS = (.03, .10, .25, .50, 1.50)


def covered_days(t, start, end, w):
    t=np.unique(t[(t>=start-w)&(t<=end+w)])
    if len(t)<4:
        return 0.
    # In the open interval (t[i], t[i+1]), the second closest exposures
    # on each side occur at t[i-1] and t[i+2]. Endpoints have zero measure.
    lo=np.maximum(np.maximum(t[1:-2],t[3:]-w),start)
    hi=np.minimum(np.minimum(t[2:-1],t[:-3]+w),end)
    return float(np.maximum(hi-lo,0.).sum())


def one(job):
    name=job[0]
    year=2000+int(name[4:6]) if name.startswith('Gaia') else int(name.split('-')[1])
    epoch=np.datetime64('2000-01-01')
    def jd(s):
        return 2451544.5+int((np.datetime64(f'{year}-{s}')-epoch).astype(int))
    start,end=jd('03-01'),jd('11-01')
    try:
        t,f,e,g,_=read_curve(job)
        t,*_=clean(t,f,e,g)
        row={'event':name,'season_d':end-start,'distinct_epochs_in_season':len(np.unique(t[(t>=start)&(t<end)]))}
        for w in WIDTHS:
            row[f'days_{w:g}']=round(covered_days(t,start,end,w),5)
        return row
    except (OSError,ValueError,KeyError,IndexError) as ex:
        return {'event':name,'error':str(ex)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,default=Path.cwd())
    p.add_argument('--out',type=Path,required=True);p.add_argument('--workers',type=int,default=3)
    a=p.parse_args()
    with mp.Pool(a.workers) as pool, a.out.open('w') as out:
        cols=['event','season_d','distinct_epochs_in_season']+[f'days_{w:g}' for w in WIDTHS]+['error']
        writer=csv.DictWriter(out,cols);writer.writeheader()
        for row in pool.imap(one,inputs(a.cache),chunksize=5):
            writer.writerow(row)
