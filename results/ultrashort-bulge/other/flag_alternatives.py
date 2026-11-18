"""Profile the strongest flags from flag_refits.py with physical binary alternatives.

Selection is defined in candidate_ids. Every eligible case is exported even if a
photometry stream fails to load. Compares all lens families on a common window.
"""
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from fetch import write_csv
from flag_refits import HERE,SEED,WORK,inputs
from peak_models import solve
from consolidate import write_gz

CACHE=Path('staging/flag-alternatives')


def candidate_ids():
    short=[];large=[]
    for p in WORK.glob('*.json'):
        r=json.loads(p.read_text())
        if r.get('status')!='ok':continue
        if r['fs_te_or_rise']<.5 and r['fs_over_both_huber']>=15:short.append(r['event'])
        if min(r['n_peak_before_3sigma'],r['n_peak_after_3sigma'])<2:continue
        if r['rho']>1 and r['fs_over_both_huber']>=25:large.append(r['event'])
    return set(short)|set(large)


def detailed(job):
    n=job[0];dest=CACHE/(n+'.json')
    if dest.exists():return
    try:
        s=json.loads((SEED/(n+'.json')).read_text())
        r=json.loads((WORK/(n+'.json')).read_text())
        # Keep the same origin convention as peak_models.solve and preserve grid
        # binary geometry. Start the PS/FS optimizations at their refined peaks.
        s.update(ps_t0=r['ps_t0'],ps_te=r['ps_te_or_rise'],ps_u0=r['ps_u0_or_fall'],
                 fs_t0=r['fs_t0'],fs_te=r['fs_te_or_rise'],fs_u0=r['fs_u0_or_fall'],fs_rho=r['rho'])
        half=float(np.clip(max(4,2.7*s['fs_te']*max(1,s['fs_rho'],s['fs_u0'])),4,80))
        bin_d=float(np.clip(min(s['fs_te'],s['ps_te'])/12,.008,.24))
        res=solve(job,s,1.15,12,True,half,bin_d)
        m=res['models']
        row=dict(event=n,status='ok',n_peak_bins=res['n_peak_window'],width_d=round(half,4),bin_d=round(bin_d,5),
                 fs_t0=round(m['fs']['t0'],6),fs_te=round(m['fs']['te_or_rise_d'],5),fs_rho=round(m['fs']['rho'],5),
                 ps_chi2=m['ps']['chi2'],fs_chi2=m['fs']['chi2'],flare_chi2=m['exp']['chi2'],binary_chi2=m['bin']['chi2'],
                 binary_te=round(m['bin']['te_or_rise_d'],5),binary_s=round(m['bin']['s'],5),binary_q=round(m['bin']['q'],6),
                 binary_angle=round(m['bin']['alpha'],5),binary_converged=m['bin']['converged'],
                 chi2_rho_le1=m['profile_small']['chi2'],chi2_rho_ge1=m['profile_large']['chi2'],
                 chi2_te_le0p5=m['profile_short']['chi2'],chi2_te_ge0p5=m['profile_long']['chi2'])
    except Exception as ex:
        row=dict(event=n,status='failed',error=repr(ex))
    tmp=dest.with_suffix('.tmp');tmp.write_text(json.dumps(row,allow_nan=False));tmp.replace(dest)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--publish-only',action='store_true')
    a=p.parse_args();CACHE.mkdir(parents=True,exist_ok=True);ids=candidate_ids()
    all_jobs={j[0]:j for j in inputs(Path.cwd()) if j[0] in ids}
    if not a.publish_only:
        priority=sorted(ids,key=lambda n:(n!='OGLE-2026-BLG-0670','2026' not in n,n))
        with ProcessPoolExecutor(a.workers) as workers:
            list(workers.map(detailed,[all_jobs[n] for n in priority if n in all_jobs]))
    rows=[json.loads(p.read_text()) for p in sorted(CACHE.glob('*.json'))]
    if rows:
        write_csv(HERE/'flag_alternatives.csv',rows)
    print(len(rows),'available out of',len(ids),flush=True)
