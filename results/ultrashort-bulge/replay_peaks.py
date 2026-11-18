"""Recompute each published peak-model chi-square from its parameter vector and raw light curve."""
import argparse
import csv
import gzip
import json
from pathlib import Path

from peak_models import solve
from run_grid import HERE, inputs
from replay_guard import external,fresh_for_inputs,verified,verify_receipt


def run(cache,out):
    cache=external(cache)
    jobs={j[0]:j for j in inputs(cache)}
    verify_receipt('tables/peak-comparisons.json.gz')
    models=json.loads(gzip.decompress((HERE/'tables/peak-comparisons.json.gz').read_bytes()))
    verify_receipt('tables/grid-fits.csv.gz')
    with gzip.open(HERE/'tables/grid-fits.csv.gz','rt',newline='') as f:
        frozen={r['event']:r['input_sha256'] for r in csv.DictReader(f)}
    for event in models:
        if event not in jobs:raise ValueError('missing frozen event: '+event)
        if not (event.startswith('OGLE-') or event.startswith('KMT-')):
            raise ValueError('derived NPZ has no frozen byte lineage: '+event)
        verified(jobs[event][1],frozen[event])
    out=fresh_for_inputs(out,cache)
    results=[]
    for event,experiments in sorted(models.items()):
        for name,run in experiments.items():
            if not run:continue
            receipt=dict(event=event,experiment=name)
            try:
                obtained=solve(jobs[event],{'fs_t0':run['reference_jd']},run['cap_factor'],
                               window=run.get('window_halfwidth_d',4.),bin_width=run.get('bin_width_d',0.),
                               evaluate=run['models'])
                distance={k:abs(v-run['models'][k]['chi2']) for k,v in obtained.items()}
                receipt.update(models=len(distance),max_abs_chi2=round(max(distance.values()),6),
                               largest_at=max(distance,key=distance.get))
                # Stored chi-squares have three decimal places; HJD rounding in the
                # constrained-profile epoch can have larger impact at a narrow caustic.
                receipt['status']='ok' if max(distance.values())<=.05 else 'check'
            except (ValueError,KeyError,IndexError,OSError) as err:
                receipt.update(status='error',detail=str(err))
            results.append(receipt)
    cols=list(dict.fromkeys(k for r in results for k in r))
    with (out/'peak-replay.csv').open('x',newline='') as f:
        w=csv.DictWriter(f,cols,lineterminator='\n');w.writeheader();w.writerows(results)
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    rows=run(args.cache,args.out);print('Replayed',len(rows),'experiments')
    for r in rows:
        if r['status']!='ok':print(r)
