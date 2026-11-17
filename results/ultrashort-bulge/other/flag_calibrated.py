"""Recompute all flagged MOA/PRIME light curves with a usable reference-flux ceiling.

Reference counts and zero point are the survey's own calibration. Add the median
off-window difference flux, then multiply by cap_factor for a sensitivity test.
No reference flux is inferred for rows without this independent constraint.
"""
import argparse
import csv
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from fetch import OUT,write_csv
from flag_refits import fit_models,SEED,WORK,inputs
from consolidate import write_gz

DIR=Path('staging/flag-calibrated')


def refit(args):
    job,c,factor=args
    dest=DIR/(job[0]+f'-{factor:.2f}.json')
    if not dest.exists():
        s=json.loads((SEED/(job[0]+'.json')).read_text())
        try:r=fit_models(job,s,c,factor)
        except Exception as ex:r=dict(event=job[0],status='error',error=repr(ex))
        tmp=dest.with_suffix('.tmp');tmp.write_text(json.dumps(r,allow_nan=False));tmp.replace(dest)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--sensitivity',action='store_true')
    a=p.parse_args();DIR.mkdir(parents=True,exist_ok=True)
    active={p.stem for p in WORK.glob('*.json') if p.stem.startswith(('MOA','PRIME'))}
    cals={r['event']:r for r in csv.DictReader((OUT/'page_calibrations.csv').open()) if r['event'] in active}
    jobs={j[0]:j for j in inputs(Path.cwd()) if j[0] in active}
    strong={r['event'] for r in csv.DictReader((OUT/'flag_alternatives.csv').open())}
    tasks=[(jobs[n],c,fac) for n,c in cals.items() if c['reference_status']=='positive reference' for fac in
           ((1.,1.15,1.35) if a.sensitivity and n in strong else (1.15,))]
    with ProcessPoolExecutor(a.workers) as pool:list(pool.map(refit,tasks))
    rows=[]
    for n,c in sorted(cals.items()):
        f=DIR/(n+'-1.15.json')
        if f.exists():rows.append(dict(reference_status=c['reference_status'],**json.loads(f.read_text())))
        else:rows.append(dict(event=n,reference_status=c['reference_status'],status='no positive independent reference'))
    write_gz(OUT/'calibrated-flag-refits.csv.gz',rows)
    if a.sensitivity:
        write_csv(OUT/'flag_calibration_sensitivity.csv',[json.loads(f.read_text()) for f in sorted(DIR.glob('*.json')) if f.name.rsplit('-',1)[0] in strong])
