"""Package the primary trial receipts and independently check the four zero-count limits."""
import argparse
import csv
import gzip
import io
import json
import math
from collections import defaultdict
from pathlib import Path

from scipy.optimize import brentq
from scipy.stats import beta


def write_csv(rows):
    buf = io.StringIO()
    writer = csv.DictWriter(buf, rows[0].keys(), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return buf.getvalue()


def package(source, out):
    summaries, trials = [], []
    for filename in sorted(source.glob('*-P*-R*.json')):
        summaries.append(json.loads(filename.read_text()))
        with filename.with_suffix('.csv').open() as f:
            trials.extend(csv.DictReader(f))
    out.mkdir(parents=True, exist_ok=True)
    (out/'ffi_recovery_stars.csv').write_text(write_csv(summaries))
    (out/'ffi_recovery_trials.csv.gz').write_bytes(gzip.compress(write_csv(trials).encode(), mtime=0))
    (out/'ffi_limits.json').write_bytes((source/'limits.json').read_bytes())


def check(where):
    def key(r):
        return r['target'],float(r['period_days']),float(r['radius_earth'])
    trial_groups = defaultdict(list)
    with gzip.open(where/'ffi_recovery_trials.csv.gz','rt') as stream:
        for row in csv.DictReader(stream):
            trial_groups[key(row)].append(row)
    with (where/'ffi_recovery_stars.csv').open() as stream:
        stars = list(csv.DictReader(stream))
    bounds = json.loads((where/'ffi_limits.json').read_text())
    assert len(stars) == len(trial_groups) == 12 and len(bounds) == 4
    assert {key(row) for row in stars} == set(trial_groups)
    for row in stars:
        trials = trial_groups[key(row)]
        assert len(trials) == int(row['trials']) == 100
        assert sorted(int(t['trial']) for t in trials) == list(range(100))
        assert all(int(t['recovered']) == int(bool(t['first_visit_index'])) for t in trials)
        assert sum(int(t['recovered']) for t in trials) == int(row['recovered'])
        assert all(0 <= float(t['impact_fraction']) <= 1 and float(t['duration_min']) > 0 for t in trials)
    for result in bounds:
        selected = [r for r in stars if key(r)[1:] == (result['period_days'], result['radius_earth'])]
        assert len({r['target'] for r in selected}) == len(selected) == result['targets'] == 3
        assert sum(int(r['possible_host']) for r in selected) == result['possible_hosts'] == 0
        assert sum(int(r['recovered']) for r in selected) == result['recovered']
        assert sum(int(r['trials']) for r in selected) == result['trials']
        rates = [int(r['recovered']) / int(r['trials']) for r in selected]
        gs = [float(r['geometry_probability']) for r in selected]
        assert math.isclose(sum(rates), result['effective_transiting_stars'])
        assert math.isclose(sum(x*y for x,y in zip(rates,gs)), result['geometry_weighted_stars'])
        lower = [beta.ppf(.025/len(selected), int(r['recovered']), 1+int(r['trials'])-int(r['recovered']))
                 if int(r['recovered']) else 0. for r in selected]
        for field,geometry in [('transiting_host_upper95',[1.]*len(gs)), ('intrinsic_host_upper95',gs)]:
            def tail(f):
                return math.prod(1-f*l*g for l,g in zip(lower,geometry))-.025
            if tail(1.) >= 0.:
                assert result[field] is None
            else:
                assert math.isclose(result[field], brentq(tail,0.,1.), rel_tol=1e-10)
    print('Reconciled 1,200 individual trials, 12 target/cell sensitivities and four pointwise limits.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record',type=Path)
    parser.add_argument('--from-run',type=Path,help='convert recovery.py intermediate files before auditing')
    args = parser.parse_args()
    if args.from_run:
        package(args.from_run,args.record)
    check(args.record)
