"""Calendar- and orbit-matched annual contrasts, then common-stage contrasts.

Descriptive tests for persistent departures, using each observation at most once.
The exact sign test uses separated 42-day block medians as the replicate.
"""
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import binomtest

from coverage import geometry_ok

BASE = Path(__file__).resolve().parent
FIELDS = {'wse':'wse_u','width':'width_u','slope':'slope_u','sos_consensus_q':None}
MAX_DAYS = 14
BLOCK_DAYS = 42


def pairs_for_pass(rows, require_stage):
    first = rows[rows.day < pd.Timestamp('2024-07-01')]
    second = rows[rows.day >= pd.Timestamp('2024-07-01')]
    if first.empty or second.empty:
        return []
    years_later = first.day + pd.DateOffset(years=1)
    gaps = np.abs((second.day.to_numpy()[None,:] - years_later.to_numpy()[:,None])/np.timedelta64(1,'D'))
    ok = gaps <= MAX_DAYS
    if require_stage:
        dh = np.abs(first.wse.to_numpy()[:,None] - second.wse.to_numpy()[None,:])
        ok &= dh <= .20
    # Rectangular assignment maximizes eligible cardinality, then minimizes offset.
    a,b = linear_sum_assignment(np.where(ok, gaps, 100000.))
    return [(first.iloc[i],second.iloc[j]) for i,j in zip(a,b) if ok[i,j]]


def main():
    manifest = json.loads((BASE/'data/manifest.json').read_text())
    pairs, summary = [], []
    for source in manifest['hydrocron']:
        with (BASE/'data'/source['file']).open() as f:
            rows = list(csv.DictReader(f))
        for gate in (1,2):
            clean = pd.DataFrame([r for r in rows if geometry_ok(r,gate)])
            if clean.empty:
                continue
            clean['day'] = pd.to_datetime(clean.time_str,utc=True).dt.normalize().dt.tz_localize(None)
            for k in [*FIELDS,*[v for v in FIELDS.values() if v]]:
                clean[k] = pd.to_numeric(clean[k],errors='coerce').mask(lambda x:x <= -1e10)
            for stage in (False,True):
                sample = [p for _,g in clean.groupby('pass_id') for p in pairs_for_pass(g,stage)]
                for field,err in FIELDS.items():
                    if stage and field in ('wse','sos_consensus_q'):
                        continue
                    observations = []
                    for a,b in sample:
                        if not np.isfinite([a[field],b[field]]).all():
                            continue
                        d = b[field]-a[field]
                        pair_sigma = np.hypot(a[err],b[err]) if err else np.nan
                        rec = dict(reach_id=source['reach_id'],max_quality=gate,stage_matched=stage,
                                   field=field,pass_id=a.pass_id,early=a.day.date(),late=b.day.date(),
                                   offset_days=(b.day-(a.day+pd.DateOffset(years=1))).days,
                                   early_crid=a.crid,late_crid=b.crid,
                                   early_value=a[field],late_value=b[field],delta=d,
                                   pair_formal_sigma=pair_sigma,
                                   stage_delta_m=b.wse-a.wse,
                                   block=(a.day-pd.Timestamp('2023-07-01')).days//BLOCK_DAYS)
                        observations.append(rec)
                        pairs.append(rec)
                    if not observations:
                        continue
                    sample_df = pd.DataFrame(observations)
                    blocks = sample_df.groupby('block').delta.median()
                    positive, negative = (blocks>0).sum(), (blocks<0).sum()
                    n = int(positive+negative)
                    summary.append(dict(reach_id=source['reach_id'],site=source['site'],max_quality=gate,
                                        stage_matched=stage,field=field,n_pairs=len(observations),blocks=n,
                                        positive_blocks=positive,median_delta=sample_df.delta.median(),
                                        same_crid_pairs=(sample_df.early_crid==sample_df.late_crid).sum(),
                                        raw_sign_p=binomtest(positive,n,.5).pvalue if n else 1.))
    s = pd.DataFrame(summary)
    # Six distinct endpoints per reach, including both stage-conditioned geometry endpoints.
    s['bonferroni_p'] = np.minimum(1, s.raw_sign_p*12*6)
    pd.DataFrame(pairs).to_csv(BASE/'annual_pairs.csv',index=False,float_format='%.9g')
    s.to_csv(BASE/'annual_screen.csv',index=False,float_format='%.9g')


if __name__ == '__main__':
    main()
