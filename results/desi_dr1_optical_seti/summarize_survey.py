"""Apply measured cframe transforms and line masks to injection ledgers."""
import argparse
import csv
import gzip
import io
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta

from vetting import reference_mask


def write_csv(file, rows):
    writer = csv.DictWriter(file, fieldnames=tuple(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


def grid_limit(counts, n, target, conf=.95):
    if not n:
        return [0.]*len(counts), None
    alpha = (1-conf)/len(counts)
    lower = [float(beta.ppf(alpha, k, n-k+1)) if k else 0. for k in counts]
    index = None
    for i in reversed(range(len(lower))):
        if lower[i] < target:
            break
        index = i
    return lower, index


def compact_rows(raw, meta, headers):
    epochs, trial = raw['epoch'], raw['trials']
    arm = raw['arm']
    petal = {x['expid']:x['petal_loc'] for x in meta['exposures']}
    us = np.array([x['rest_A'] for x in trial])
    coverage = np.array([x['eligible'] for x in trial])
    eligible = coverage.copy()
    fluxes = raw['flux17']
    reasons = np.asarray([t['reasons'] for t in trial],dtype=object)
    for i, e in enumerate(epochs):
        fac = headers[e['expid'],arm+str(petal[e['expid']])]
        d = e['doppler']
        fwhm = np.array([t['fwhm_A'][i] if t['eligible'] else np.nan for t in trial])
        excluded = reference_mask(us*d, fwhm, fac, d, e['rv_error'])
        eligible &= ~excluded
        for k in range(len(fluxes)):
            idx = np.flatnonzero(reasons[:,k]=='recovered')
            if len(idx):
                fit = np.array([trial[j]['fitted'][k][i] for j in idx])
                rejected = reference_mask(fit[:,0],fit[:,1],fac,d,e['rv_error'])
                reasons[idx[rejected],k] = 'fitted_reference_veto'
    reasons[coverage & ~eligible,:] = 'injected_reference_veto'
    codes = dict(recovered='R',missed_detection='M',baseline_coincidence='B',outside_joint_coverage='G',
                 injected_reference_veto='V',fitted_reference_veto='F')
    ledger = [dict(targetid=raw['targetid'],arm=arm,draw=j,rest_A=repr(float(w)),
                   states=''.join(codes[r] for r in reasons[j])) for j,w in enumerate(us)]
    counts = (reasons=='recovered').sum(0).tolist()
    n = len(trial)
    n_good = int(eligible.sum())
    row = dict(targetid=raw['targetid'],arm=arm,seed=raw['seed'],expids='|'.join(str(e['expid']) for e in epochs),
               rest_min_A=raw['range_rest_A'][0],rest_max_A=raw['range_rest_A'][1],trials=n,
               n_joint_coverage=int(coverage.sum()),n_joint_reference_eligible=n_good)
    for threshold in (.50,.80,.90,.95):
        for nc,denominator in ((n,'all'),(n_good,'eligible')):
            _, level = grid_limit(counts,nc,threshold)
            row[f'flux{round(threshold*100)}_{denominator}_17'] = '' if level is None else fluxes[level]
    low_all,_ = grid_limit(counts,n,.95)
    low_ok,_ = grid_limit(counts,n_good,.95)
    levels = [dict(targetid=raw['targetid'],arm=arm,level=k,flux17=fluxes[k],recovered=c,
                   n_all=n,n_eligible=n_good,lower95_all=low_all[k],lower95_eligible=low_ok[k])
              for k,c in enumerate(counts)]
    return row,levels,ledger


def run(scan_dir, inject_dir, header_files, outpath, tag='*'):
    meta = {}
    for path in sorted(Path(scan_dir).glob('scan-*.jsonl')):
        for line in path.read_text().splitlines():
            row=json.loads(line)['meta']
            meta[row['targetid']]=row
    headers={}
    for path in header_files:
        with Path(path).open() as f:
            for h in csv.DictReader(f):
                k = (int(h['expid']),h['camera'])
                if k in headers and float(h['heliocor'])!=headers[k]:
                    raise ValueError('conflicting HELIOCOR')
                headers[k]=float(h['heliocor'])
    rows,levels,ledger=[],[],[]
    seen=set()
    for path in sorted(Path(inject_dir).glob(f'injections-{tag}.jsonl.gz')):
        with gzip.open(path,'rt') as f:
            for line in f:
                raw=json.loads(line)
                k=(raw['targetid'],raw['arm'])
                if k in seen:
                    raise ValueError('duplicate target and arm across partitions')
                seen.add(k)
                r,l,t=compact_rows(raw,meta[k[0]],headers)
                rows.append(r)
                levels.extend(l)
                ledger.extend(t)
    if not rows:
        raise ValueError('no injection records found')
    out = Path(outpath)
    out.mkdir(exist_ok=True,parents=True)
    for name, rr in [('flux_limits.csv',rows),('recovery_grid.csv',levels)]:
        rr.sort(key=lambda r: (r['targetid'],r['arm']))
        with (out/name).open('w',newline='') as f:
            write_csv(f,rr)
    ledger.sort(key=lambda r:(r['targetid'],r['arm'],r['draw']))
    with (out/'injection_trials.csv.gz').open('wb') as f:
        with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0) as compressed:
            with io.TextIOWrapper(compressed,encoding='utf-8',newline='') as t:
                write_csv(t,ledger)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scan_dir')
    parser.add_argument('inject_dir')
    parser.add_argument('outpath')
    parser.add_argument('headers',nargs='+')
    parser.add_argument('--tag',default='*')
    args=parser.parse_args()
    run(args.scan_dir,args.inject_dir,args.headers,args.outpath,args.tag)
