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
        # Rebuild the compatible-host numerator from the original blind search.
        m = accepted[r['target']]
        period, radius = key(r)[1:]
        mass, rs = m['mass_solar'], m['radius_solar']*6.957e8
        a = (1.3271244e20*mass*(period*86400/(2*math.pi))**2)**(1/3)
        g = (rs + radius*6.3781e6)/a
        assert math.isclose(float(r['geometry_probability']), g, rel_tol=1e-12)
        max_width = period*math.asin(g)/math.pi
        cycles = math.floor((m['end_btjd']-m['start_btjd'])/period)
        tolerance = max_width + 2*m['cadence_s']/86400
        eligible = any(cycles*abs(x['period_days']-period) <= tolerance
                       for x in m['search']['periodic'])
        assert int(eligible) == int(r['possible_host'])
        key_bytes = f"{r['target']}|{period:.12g}|{radius:.12g}|337109".encode()
        seed = int.from_bytes(hashlib.sha256(key_bytes).digest()[:8], 'big')
        rng = np.random.default_rng(seed)
        for trial in sorted(group, key=lambda t: int(t['trial'])):
            epoch = m['start_btjd']+rng.random()*period
            impact = rng.random()
            assert abs(epoch-float(trial['epoch_btjd'])) < 1e-10
            assert impact == float(trial['impact_fraction'])
            fraction = g*math.sqrt((1-impact**2)/(1-(g*impact)**2))
            duration = period*1440*math.asin(fraction)/math.pi
            assert math.isclose(duration, float(trial['duration_min']), rel_tol=1e-10)
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
    # Reconcile the separate aggregation method using only the trial ledger.
    for cell in json.loads((root/'weighted_bounds.json').read_text()):
        groups = {target: group for (target, period, radius), group in trials.items()
                  if (period, radius) == (cell['period'], cell['radius_earth'])}
        assert set(groups) == set(accepted)
        assert cell['targets'] == len(groups) and cell['possible'] == 0
        summary = {r['target']: r for r in summaries if key(r)[1:] == (cell['period'], cell['radius_earth'])}
        assert all(int(r['possible_host']) == 0 for r in summary.values())
        for label in ('transiting', 'intrinsic'):
            mean, variance_scale = 0., 0.
            for target, group in groups.items():
                weight = 1. if label == 'transiting' else float(summary[target]['geometry_probability'])
                mean += weight*sum(int(t['recovered']) for t in group)/len(group)
                variance_scale += weight**2/len(group)
            penalty = math.sqrt(-.5*variance_scale*math.log(.005))
            floor = max(0., mean-penalty)
            result = cell[label]
            assert math.isclose(mean, result['effective_stars'])
            assert math.isclose(penalty, result['mc_penalty'])
            assert math.isclose(floor, result['mc_lower'])
            f = result['upper95']
            tail = (1.-(1. if f is None else f)*floor/len(groups))**len(groups)
            assert tail >= .045 if f is None else math.isclose(tail, .045, rel_tol=1e-10)
    print(f'Checked {len(roster)} parent records, {len(accepted)} independent stars and {sum(map(len,trials.values()))} injections.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record',type=Path)
    check(parser.parse_args().record)
