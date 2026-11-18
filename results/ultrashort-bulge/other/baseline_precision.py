"""Measure the precision and temporal scatter of the off-event flux calibration.

The bootstrap resamples populated nights. It measures sampling precision only.
It cannot measure an absolute zero-point error or independent source blending.
"""
import csv
import json
from pathlib import Path

import numpy as np

from fetch import OUT, write_csv
from flag_refits import SEED, inputs, window


def summarize(job, calibration):
    seed = json.loads((SEED / (job[0] + '.json')).read_text())
    _, meta, (time, flux, _, _) = window(job, seed)
    origin, half, _, _, num_off, _, raw_median = meta
    center = seed['ps_t0'] - origin
    start, stop = min(0., center) - half, max(0., center) + half
    off = (time < start - half / 2) | (time > stop + half / 2)
    if off.sum() < 3:
        off = (time < start) | (time > stop)
    assert off.sum() == num_off
    row = dict(event=job[0], reference_status=calibration['reference_status'],
               off_epochs=int(num_off), off_epoch_median_counts=round(raw_median, 3))
    if not num_off:
        return row
    time, flux = time[off] + origin, flux[off]
    nights = np.floor(time).astype(int)
    split = np.flatnonzero(np.diff(nights)) + 1
    medians = np.array([np.median(chunk) for chunk in np.split(flux, split)])
    row['off_nights'] = len(medians)
    row['off_nightly_median_counts'] = round(float(np.median(medians)), 3)
    q16, q84 = np.quantile(medians, [.16, .84])
    row['off_nightly_p16_counts'] = round(q16, 3)
    row['off_nightly_p84_counts'] = round(q84, 3)
    if calibration['reference_status'] != 'positive reference':
        return row
    ref = float(calibration['additive_reference_counts'])
    baseline = ref + raw_median
    row['source_baseline_counts'] = round(baseline, 3)
    if len(medians) < 3 or baseline <= 0:
        return row
    rng = np.random.default_rng(614)
    samples = [np.median(rng.choice(medians, len(medians), replace=True)) for _ in range(400)]
    lo, hi = np.quantile(samples, [.16, .84])
    row['nightly_bootstrap_p16_counts'] = round(lo + ref, 3)
    row['nightly_bootstrap_p84_counts'] = round(hi + ref, 3)
    row['p84_minus_p16_over_twice_baseline'] = round((hi - lo) / (2 * baseline), 6)
    row['nightly_p84_minus_p16_over_twice_baseline'] = round((q84 - q16) / (2 * baseline), 6)
    return row


if __name__ == '__main__':
    directory = Path('staging/flag-refits')
    calibrations = {r['event']: r for r in csv.DictReader((OUT / 'page_calibrations.csv').open())
                    if (directory / (r['event'] + '.json')).exists()}
    jobs = sorted((j for j in inputs(Path.cwd()) if j[0] in calibrations), key=lambda j: j[0])
    write_csv(OUT / 'baseline_precision.csv', [summarize(j, calibrations[j[0]]) for j in jobs])
