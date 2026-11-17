"""Expand binary and flare fits for the best short or large-source hypotheses.

Reuse the peak-window analysis' exact reference epoch, binning and stream caps.
All reported differences compare the same measurements and error weights.
"""
import argparse
import json
import multiprocessing as mp
from pathlib import Path

from peak_models import solve
from run_grid import inputs

HERE = Path(__file__).resolve().parent


def priority(record):
    m = record['models']
    if 'profile_small' not in m:
        return False
    single = min(v['chi2'] for k, v in m.items() if k == 'fs' or k.startswith('profile_'))
    small_cost = m['profile_small']['chi2'] - single
    long_cost = m['profile_long']['chi2'] - single
    return ((small_cost > 7 or long_cost > 2) and
            m['bin']['chi2'] > single - 20 and m['exp']['chi2'] > single - 8)


def run(args):
    job, seed, record, out = args
    seed = dict(seed, fs_t0=record['reference_jd'])
    new = solve(job, seed, record['cap_factor'], 8, False,
                record.get('window_halfwidth_d', 4.), record.get('bin_width_d', 0.),
                previous=record['models'])
    assert new['n_peak_window'] == record['n_peak_window']
    p = out / (job[0] + '.json')
    tmp = p.with_suffix('.tmp')
    tmp.write_text(json.dumps(new, allow_nan=False) + '\n')
    tmp.replace(p)
    return job[0]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path.cwd())
    parser.add_argument('--workers', type=int, default=2)
    a = parser.parse_args()
    out = a.cache / 'staging/expanded'
    out.mkdir(exist_ok=True, parents=True)
    jobs = {j[0]: j for j in inputs(a.cache)}
    selected = []
    for folder in (a.cache/'staging/comparisons', HERE/'peaks'):
        for path in folder.glob('*.json'):
            r = json.loads(path.read_text())
            if 'models' not in r or (out/path.name).exists() or not priority(r) or path.stem not in jobs:
                continue
            bound = a.cache/'staging/bound'/path.name
            seed = json.loads((bound if bound.is_file() else a.cache/'staging/screen'/path.name).read_text())
            selected.append((jobs[path.stem], seed, r, out))
    print('Queued', len(selected), flush=True)
    with mp.Pool(a.workers) as pool:
        for event in pool.imap_unordered(run, selected):
            print(event, flush=True)
