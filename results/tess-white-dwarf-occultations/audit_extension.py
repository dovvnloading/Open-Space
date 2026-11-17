"""Audit extension membership, every injection, archive hashes and Poisson-binomial upper limits."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import beta


def table(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt') as stream:
        return list(csv.DictReader(stream))


def check(root):
    selection=json.loads((root/'extension_selection.json').read_text())
    parent=root/'extension_parent.csv'
    assert hashlib.sha256(parent.read_bytes()).hexdigest()==selection['roster_sha256']
    roster=table(parent)
    assert Counter(r['selection'] for r in roster)==selection['decisions']
    selected={r['GaiaEDR3'] for r in roster if r['selection']=='request latest visit'}
    outcomes=json.loads((root/'extension_measurements.json').read_text())
    assert len(selected)==len(outcomes) and selected=={r['source_id'] for r in outcomes}
    accepted={r['source_id']:r for r in outcomes if r['status']=='included'}
    receipts={r['file']:r for r in json.loads((root/'extension_receipts.json').read_text())}
    with gzip.open(root/'extension_catalogue_responses.json.gz','rt') as stream:
        catalogue_inputs=json.load(stream)
    for filename,payload in catalogue_inputs.items():
        assert hashlib.sha256(payload.encode('utf-8')).hexdigest()==receipts[filename+'.json']['sha256']
    assert {r['target'] for r in table(root/'extension_visits.csv')}==set(accepted)
    for r in accepted.values():
        token=f"{r['source_id']}-s{r['sector']}.fits.json"
        assert receipts[token]['sha256']==r['fits_sha256'] and receipts[token]['http_status']==200
        assert r['min_target_fraction']>=.8 and r['gaia_without_G']==0 and r['quality0_points']>=100
        assert abs(r['cadence_s']-200)<1 and 56<=r['sector']<=108
    def key(r):
        return r['target'],float(r['period_days']),float(r['radius_earth'])
    trials=defaultdict(list)
    for r in table(root/'extension_recovery_trials.csv.gz'):
        trials[key(r)].append(r)
    summaries=table(root/'extension_recovery_stars.csv')
    assert len(trials)==len(summaries)==4*len(accepted)
    assert {key(r) for r in summaries}==set(trials)
    for r in summaries:
        group=trials[key(r)]
        assert sorted(int(t['trial']) for t in group)==list(range(50))
        assert sum(int(t['recovered']) for t in group)==int(r['recovered'])
        assert int(r['trials'])==50 and r['target'] in accepted
        assert all(int(t['recovered'])==int(bool(t['first_visit_index'])) for t in group)
    for cell in json.loads((root/'extension_limits.json').read_text()):
        rows=[r for r in summaries if key(r)[1:]==(cell['period_days'],cell['radius_earth'])]
        n=len(rows)
        assert n==cell['targets']==len(accepted)
        count=sum(int(r['possible_host']) for r in rows)
        assert count==cell['possible_hosts']
        assert sum(int(r['recovered']) for r in rows)==cell['recovered']
        q=np.array([int(r['recovered'])/50 for r in rows])
        g=np.array([float(r['geometry_probability']) for r in rows])
        assert math.isclose(q.sum(),cell['effective_transiting_stars'])
        assert math.isclose((q*g).sum(),cell['geometry_weighted_stars'])
        lo=np.array([beta.ppf(.025/n,k,51-k) if k else 0. for k in (int(r['recovered']) for r in rows)])
        for prop,field in ((lo,'transiting_host_upper95'),(lo*g,'intrinsic_host_upper95')):
            f=cell[field]
            p=prop*(1. if f is None else f)
            # EGF for odds gives exactly P(K=j) / P(K=0) as coefficient j.
            odds=p/(1.-p)
            coef=np.array([1.]+[0.]*count)
            for odd in odds:
                coef[1:]+=odd*coef[:-1]
            tail=math.prod(1.-p)*coef.sum()
            if f is None:
                assert tail >= .025-1e-12
            else:
                assert math.isclose(tail,.025,rel_tol=1e-8)
    print(f'Checked {len(roster)} parent records, {len(accepted)} independent stars and {sum(map(len,trials.values()))} injections.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record',type=Path)
    check(parser.parse_args().record)
