"""Profile primary/near-threshold leads under survey-calibrated source flux bounds."""
import csv
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from fetch import OUT,write_csv
from flag_refits import inputs,SEED
from peak_models import solve

DEST=Path('staging/flag-calibrated-binary')


def work(task):
    job,r=task
    n=job[0];dest=DEST/(n+'.json')
    if dest.exists():return
    s=json.loads((SEED/(n+'.json')).read_text())
    s.update(ps_t0=r['ps_t0'],ps_te=r['ps_te_or_rise'],ps_u0=r['ps_u0_or_fall'],
             fs_t0=r['fs_t0'],fs_te=r['fs_te_or_rise'],fs_u0=r['fs_u0_or_fall'],fs_rho=r['rho'])
    hw=float(np.clip(max(4,2.8*s['fs_te']*max(1,s['fs_rho'],s['fs_u0'])),4,70))
    d=float(np.clip(min(s['fs_te'],s['ps_te'])/12,.008,.20));cap=r['cap_factor']*r['baseline_counts']
    result=solve(job,s,r['cap_factor'],32,True,hw,d,source_flux_caps=[cap]);m=result['models']
    row=dict(event=n,cap_counts=round(cap,3),quiet_baseline_counts=r['baseline_counts'],bin_d=round(d,5),
             window_halfwidth=round(hw,4),n_bins=result['n_peak_window'],binary_starts=32,
             fs_t0=round(m['fs']['t0'],6),fs_te=round(m['fs']['te_or_rise_d'],6),fs_u0=round(m['fs']['u0_or_fall_d'],6),
             fs_rho=round(m['fs']['rho'],6),source_fraction=round(m['fs']['fs'][0]/r['baseline_counts'],6),
             ps_chi2=m['ps']['chi2'],fs_chi2=m['fs']['chi2'],exp_chi2=m['exp']['chi2'],
             binary_chi2=m['bin']['chi2'],binary_q=round(m['bin']['q'],6),binary_s=round(m['bin']['s'],6),
             binary_converged=m['bin']['converged'],chi2_te_le_half=m['profile_short']['chi2'],
             chi2_te_ge_half=m['profile_long']['chi2'],chi2_rho_le_one=m['profile_small']['chi2'])
    dest.write_text(json.dumps(row)+'\n')


if __name__=='__main__':
    DEST.mkdir(parents=True,exist_ok=True)
    cals=[json.loads(p.read_text()) for p in Path('staging/flag-calibrated').glob('*-1.15.json')]
    # Add all high-score calibrated large-source solutions with two points on
    # both slopes, the best calibrated hours case, and the known 0298 cross-check.
    rows=[r for r in cals if r.get('status')=='ok' and (r['event']=='MOA-2025-BLG-0298'
          or (r['rho']>1 and r['fs_over_both_huber']>13 and min(r['n_peak_before_3sigma'],r['n_peak_after_3sigma'])>=2))]
    jobs={j[0]:j for j in inputs(Path.cwd())}
    with ProcessPoolExecutor(4) as pool:list(pool.map(work,[(jobs[r['event']],r) for r in rows]))
    write_csv(OUT/'calibrated_alternatives.csv',[json.loads(p.read_text()) for p in sorted(DEST.glob('*.json'))])
