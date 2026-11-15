"""Reproduce the closest-approach statistics from the archived state vectors."""
from datetime import datetime, timedelta
from pathlib import Path
import csv
import json
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq

DATA = Path(__file__).resolve().parent / 'data'
AU = 149597870.7
DAY = 86400.0
MJD0 = datetime(1858, 11, 17)


def approach(t, asteroid, comet):
    # The elapsed-time coordinate avoids numerical cancellation in root finding.
    elapsed = t - t[0]
    rel = asteroid - comet
    spline = CubicHermiteSpline(elapsed, rel[:, :3], rel[:, 3:])
    radial = lambda d: np.dot(spline(d), spline(d, 1))
    roots = [brentq(radial, a, b) for a, b in zip(elapsed, elapsed[1:]) if radial(a) < 0 < radial(b)]
    if len(roots) != 1:
        raise ValueError(f'Expected one interior minimum, found {len(roots)}')
    u = roots[0]
    speed = np.linalg.norm(spline(u, 1)) * AU / DAY
    miss = np.linalg.norm(spline(u)) * AU
    return {'time_tdb': (MJD0 + timedelta(days=t[0] + u)).isoformat(timespec='milliseconds'),
            'distance_km': round(miss, 2), 'speed_km_s': round(speed, 6)}


def main():
    table = np.loadtxt(DATA / 'jpl-vectors.csv', delimiter=',', skiprows=1)
    independent = np.loadtxt(DATA / 'esa-propagation.csv', delimiter=',', skiprows=1)
    if not np.allclose(table[:, 0], independent[:, 0], atol=1e-8, rtol=0):
        raise ValueError('Incompatible ephemeris times')
    result = {'jpl': approach(table[:, 0], table[:, 7:], table[:, 1:7]),
              'esa_asteroid_and_jpl_comet': approach(table[:, 0], independent[:, 1:], table[:, 1:7])}
    with open(DATA / 'jpl-observer.csv') as f:
        selected = [r for r in csv.DictReader(f) if r['utc'] == '2019-Oct-28 03:00']
    uncertainty = {}
    for r in selected:
        sky = np.deg2rad(float(r['sky_rss_3sigma_arcsec']) / 3600) * float(r['delta_au']) * AU
        rng = float(r['range_3sigma_km'])
        uncertainty[r['body']] = {'sky_rss_km': round(sky, 2), 'range_km': rng,
                                 'quadrature_km': round(float(np.hypot(sky, rng)), 2)}
    result['hourly_formal_three_sigma_scales'] = uncertainty
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
