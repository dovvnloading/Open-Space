"""Measure single-event recovery for specified dust-like profiles in the S110 visit."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.stats import beta

from wd_search import flatten, search

SEED = 11055997020
N_EPOCHS = 200
DT = 200 / 86400.
DEPTHS = (.10, .50)
DURATIONS = (4., 16.)


def trapezoid_integral(x):
    """Antiderivative of a unit-height, unit-duration symmetric trapezoid.

    x is time from first contact in units of full duration. Each linear
    ramp takes one quarter of the duration.
    """
    rise = np.clip(x, 0, .25)
    flat = np.clip(x-.25, 0, .50)
    fall = np.clip(x-.75, 0, .25)
    return 2*rise**2 + flat + fall - 2*fall**2


def cell(args):
    t, f, baseline, dilution, epochs, raw_hits, depth, minutes = args
    duration = minutes / 1440.
    trials = []

    def matched(hits, center):
        # An above-threshold isolated box must cover the injected midpoint
        # within half the injected duration or one exposure for timestamp error.
        return any(abs(r['epoch_btjd']-center) <= max(duration/2, DT) for r in hits)

    for i, center in enumerate(epochs):
        mid = (t-center)/duration + .5
        h = DT/(2*duration)
        loss = duration/DT * (trapezoid_integral(mid+h)-trapezoid_integral(mid-h))
        live = np.max(loss) > 0
        preexisting = matched(raw_hits, center)
        found = False
        if live and not preexisting:
            curve = f - depth*baseline*dilution*loss
            found = matched(search(t, curve, 200)['isolated'], center)
        trials.append(dict(trial=i, epoch_btjd=center, depth=depth, duration_min=minutes,
                           observed=int(live), preexisting_match=int(preexisting), recovered=int(found)))
    return trials


def main(series, out, workers):
    with gzip.open(series, 'rt') as f:
        rows = [r for r in map(json.loads, f) if r.get('min_target_fraction',0) >= .80]
    t, f, dilution = np.array([[r[k] for r in rows] for k in ('time','flux','target_fraction')])
    raw = search(t, f, 200)['isolated']
    base = flatten(t, f)
    base[~np.isfinite(base)] = np.nanmedian(base)
    rng = np.random.default_rng(SEED)
    epochs = rng.uniform(t.min()-DT/2, t.max()+DT/2, size=N_EPOCHS)
    arguments = [(t, f, base, dilution, epochs, raw, depth, width) for depth in DEPTHS for width in DURATIONS]
    with ProcessPoolExecutor(workers) as pool:
        cells = list(pool.map(cell, arguments))
    table = io.StringIO(newline='')
    writer = csv.DictWriter(table, fieldnames=list(cells[0][0]), lineterminator='\n')
    writer.writeheader()
    summary = []
    for trials in cells:
        writer.writerows(trials)
        n = len(trials)
        k = sum(r['recovered'] for r in trials)
        summary.append(dict(duration_min=trials[0]['duration_min'], depth=trials[0]['depth'], trials=n,
                            observed=sum(r['observed'] for r in trials), recovered=k,
                            preexisting_matches=sum(r['preexisting_match'] for r in trials),
                            mc_interval_95=[float(beta.ppf(.025,k,n-k+1)) if k else 0.,
                                            float(beta.ppf(.975,k+1,n-k)) if k<n else 1.]))
    out.mkdir(exist_ok=True, parents=True)
    trials_path = out/'single_dip_trials.csv.gz'
    trials_path.write_bytes(gzip.compress(table.getvalue().encode(),mtime=0))
    report = dict(seed=SEED, epochs_per_cell=N_EPOCHS, first_epoch_allowed=t.min()-DT/2,
                  last_epoch_allowed=t.max()+DT/2,
                  frames_sha256=hashlib.sha256(series.read_bytes()).hexdigest(),
                  trials_sha256=hashlib.sha256(trials_path.read_bytes()).hexdigest(), cells=summary)
    (out/'single_dip_recovery.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('series', type=Path)
    p.add_argument('out', type=Path)
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args()
    main(args.series, args.out, args.workers)
