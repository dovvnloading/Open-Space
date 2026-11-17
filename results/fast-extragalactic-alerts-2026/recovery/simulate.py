"""Measure conditional pair-selection recovery on the actual alert epochs."""
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
AB_ZERO_NJY = 31.4
MAG_FACTOR = 2.5 / math.log(10)
TRIALS = 16

def load(survey):
    pattern = 'rubin/data/day_*.csv.gz' if survey=='rubin' else 'ztf/data/ztf/*_detections.csv.gz'
    visit_data = defaultdict(dict)
    first = {}
    for path in sorted(ROOT.glob(pattern)):
        with gzip.open(path,'rt') as stream:
            for row in csv.DictReader(stream):
                oid, t = row['oid'], float(row['mjd'])
                first[oid] = min(t, first.get(oid,t))
                if survey=='rubin':
                    band,visit = row['band'],row['visit']
                    sigma = float(row['psffluxerr'] or 'nan')
                else:
                    band,visit = row['fid'],row['mjd']
                    depth = float(row['diffmaglim'] or 'nan')
                    sigma = 10**(0.4*(AB_ZERO_NJY-depth))/5 if 0<depth<40 else math.nan
                if sigma>0 and math.isfinite(sigma):
                    coord=(band,visit)
                    old=visit_data[oid].get(coord)
                    if old is None or sigma>old[1]:
                        visit_data[oid][coord]=(t,sigma)
    if not first:raise ValueError('Missing data for '+survey)
    return first,visit_data

def run(survey,seed):
    first,visit_data=load(survey)
    ids=sorted(first)
    groups=[]; phase=[]; errors=[]; p1=[]; p2=[]
    for idx,oid in enumerate(ids):
        by_band=defaultdict(list)
        for (band,visit),(t,sigma) in sorted(visit_data[oid].items()):
            index=len(phase); phase.append(t-first[oid]); errors.append(sigma); groups.append(idx)
            by_band[band].append((t,index))
        for visits in by_band.values():
            visits.sort()
            for i,(start,k) in enumerate(visits):
                for end,l in visits[i+1:]:
                    if end-start>1:break
                    if end-start>=0.02:p1.append(k);p2.append(l)
    group,phase,sigma=map(np.array,(groups,phase,errors))
    left,right=map(lambda x:np.array(x,dtype=int),(p1,p2))
    gap=phase[right]-phase[left]
    capable=len(np.unique(group[left]))
    results=[]
    for peak in (18,20,22):
        for rate in (0,2,4):
            rng=np.random.default_rng(seed)
            recovered=0
            for _ in range(TRIALS):
                peak_time=rng.uniform(-0.5,0.5,len(ids))
                mag=peak+rate*abs(phase-peak_time[group])
                flux=10**(0.4*(AB_ZERO_NJY-mag))+rng.normal(0,sigma)
                ok=(flux[left]>=5*sigma[left])&(flux[right]>=5*sigma[right])
                a,b=left[ok],right[ok]
                dm=2.5*np.abs(np.log10(flux[a]/flux[b]))
                dmerr=MAG_FACTOR*np.hypot(sigma[a]/flux[a],sigma[b]/flux[b])
                keep=(dm>=0.3)&(dm>gap[ok])&(dm>=3*dmerr)
                recovered+=len(np.unique(group[a[keep]]))
            results.append(dict(survey=survey,peak_AB=peak,slope_mag_day=rate,
                templates=len(ids),timing_capable=capable,injections=TRIALS*len(ids),recovered=recovered,
                percent_all=100*recovered/(TRIALS*len(ids)),percent_timing_capable=100*recovered/(TRIALS*capable)))
    return results

if __name__=='__main__':
    rows=run('rubin',20261117)+run('ztf',20261118)
    with (ROOT/'recovery/results.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    print(json.dumps(rows,indent=2))
