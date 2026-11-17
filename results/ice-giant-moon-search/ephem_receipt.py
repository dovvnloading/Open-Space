"""Freeze sampled spacecraft observables and SHA-256 receipts for the consumed responses."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def freeze(ephemeris, cache, prefix, kernel=None):
    with np.load(ephemeris) as f:
        ids, mjds, obs=f['ids'],f['mjd'],f['values']
    with open(prefix+'_observables.csv','w',newline='') as stream:
        out=csv.writer(stream)
        out.writerow(['body','mjd_utc','ra_icrf_deg','dec_icrf_deg','diam_arcsec','r_au','delta_au'])
        for body,alltimes in zip(ids,obs):
            for date,row in zip(mjds,alltimes):
                out.writerow([body,f'{date:.10f}',*[f'{x:.12g}' if np.isfinite(x) else '' for x in row]])
    # Match the response's actual JD list, avoiding unrelated cache entries for
    # the same target. These short pilot grids fit in one Horizons request each.
    with open(prefix+'_sources.csv','w',newline='') as stream:
        out=csv.writer(stream);out.writerow(['body','returned_identity','source_key','sha256'])
        for body in ids:
            if body==75052 and kernel is not None:
                out.writerow([body,'URA184 segment 7 -> 75052 and 7 -> 799',Path(kernel).name,digest(kernel)])
                continue
            hits=[]
            for path in Path(cache).glob(f'{body}-*.json'):
                raw=json.loads(path.read_bytes())['result']
                rows=raw.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines()
                dates=np.array([float(r.split(',')[0])-2400000.5 for r in rows])
                for origin in range(0,len(mjds),60):
                    chunk=mjds[origin:origin+60]
                    if len(rows)==len(chunk) and np.allclose(dates,chunk,rtol=0,atol=1e-8):
                        identity=raw.split('Target body name:')[1].splitlines()[0].strip()
                        hits.append((origin,[body,identity,path.stem,digest(path)]))
            if sorted(i for i,_ in hits)!=list(range(0,len(mjds),60)):
                raise RuntimeError(f'incomplete or ambiguous response set for {body}')
            out.writerows(row for _,row in sorted(hits))



if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ['ephemeris','cache','prefix']:
        p.add_argument(arg)
    p.add_argument('--kernel')
    a=p.parse_args();freeze(a.ephemeris,a.cache,a.prefix,a.kernel)
