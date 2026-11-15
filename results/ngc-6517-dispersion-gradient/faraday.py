"""Check a common RM-versus-DM relation for the four polarized pulsars.

Place alongside fast.csv from this result. Raw and corrected RMs in rad m^-2
share one measurement; these are sensitivity checks on the same data.
"""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parent


def load(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def fit(dm, rm, error):
    origin = np.average(dm, weights=error**-2)
    design = np.column_stack([np.ones(len(dm)), dm - origin])
    coefficient = np.linalg.lstsq(design / error[:, None], rm / error, rcond=None)[0]
    residual = rm - design @ coefficient
    statistic = np.sum((residual / error)**2)
    dof = len(dm) - len(coefficient)
    return dict(slope=float(coefficient[1]), dm_origin=origin, rm_at_origin=coefficient[0],
                residuals=residual.tolist(), chi2=statistic, dof=dof, formal_p=chi2.sf(statistic, dof))


def measure():
    rows = load(ROOT / 'faraday.csv')
    labels = [row['pulsar'] for row in rows]
    fast = {row['pulsar']: row for row in load(ROOT / 'fast.csv')}
    dm = np.array([float(fast[label]['dm_yin']) for label in labels])
    dm_err = np.array([float(fast[label]['dm_yin_err']) for label in labels])
    corrected, raw, error = [np.array([float(r[key]) for r in rows])
                             for key in ('rm_corrected', 'rm_raw', 'rm_err')]
    keep = np.arange(4) != 2
    omitted = fit(dm[keep], corrected[keep], error[keep])
    # Unknown extra variance is a sensitivity parameter, not an estimated error.
    def disagreement(floor):
        return fit(dm, corrected, np.hypot(error, floor))['formal_p'] - 0.05
    pairs = {}
    for i, j in ((0, 3), (1, 3), (2, 1)):
        delta_rm = corrected[i] - corrected[j]
        delta_dm = dm[i] - dm[j]
        pairs[labels[i][-1] + '-' + labels[j][-1]] = dict(
            delta_rm=delta_rm, rm_error_independent=np.hypot(error[i], error[j]),
            delta_dm=delta_dm, dm_error_independent=np.hypot(dm_err[i], dm_err[j]))
    return dict(pulsars=labels, corrected=fit(dm, corrected, error), raw=fit(dm, raw, error),
                omit_C=omitted, independent_rm_floor_for_p_0_05=brentq(disagreement, 0, 50),
                pairs=pairs)


if __name__ == '__main__':
    (ROOT / 'faraday.json').write_text(json.dumps(measure(), indent=2, allow_nan=False) + '\n')
