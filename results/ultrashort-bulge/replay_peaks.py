"""Recompute each published peak-model chi-square from its parameter vector and raw light curve."""
import argparse
import csv
import gzip
import json
from pathlib import Path

from peak_models import solve
from run_grid import HERE, inputs


def run(cache):
    jobs={j[0]:j for j in inputs(cache)}
    models=json.loads(gzip.decompress((HERE/'tables/peak-comparisons.json.gz').read_bytes()))
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
    with (HERE/'tables/peak-replay.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,cols,lineterminator='\n');w.writeheader();w.writerows(results)
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,default=Path.cwd());args=p.parse_args()
    rows=run(args.cache);print('Replayed',len(rows),'experiments')
    for r in rows:
        if r['status']!='ok':print(r)
