"""Recalculate ballistic miss distances with measured local proton flow vectors."""
import csv
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
AU_KM = 149597870.7
WINDOW_HOURS = (1/6, 1, 3, 6, 12, 24, 48, 72)
# Independent 400 km/s radial screening window centers, before examining plasma.
EVENTS = {'A': (2456690.16834153, 67.184), 'B': (2455149.67982216, 66.184)}


def ephemeris(probe, body):
    table = np.loadtxt(ROOT / f'ephemeris_{probe}_{body}.csv', delimiter=',', skiprows=1)
    return CubicHermiteSpline(table[:, 0], table[:, 1:4], table[:, 4:7], extrapolate=False)


def read_wind(probe, tdb_minus_utc):
    rows = list(csv.DictReader((ROOT / f'wind_{probe}.csv').open()))
    times = np.asarray([datetime.fromisoformat(r['utc']).timestamp() for r in rows])
    jd = (times + tdb_minus_utc) / 86400 + 2440587.5
    values = np.asarray([[float(row[k]) if row[k] else np.nan for k in
                          ['Vr_kms','Vt_kms','Vn_kms']] for row in rows])
    return times, jd, values


def miss(craft, comet, receipt_jd, velocity):
    position = craft(receipt_jd)
    radial = position / np.linalg.norm(position)
    tangent = np.cross([0, 0, 1], radial)
    tangent /= np.linalg.norm(tangent)
    normal = np.cross(radial, tangent)
    velocity_xyz = velocity @ np.array([radial, tangent, normal]) * 86400 / AU_KM

    def separation(lag):
        return np.linalg.norm(position - comet(receipt_jd - lag) - lag * velocity_xyz)

    optimum = minimize_scalar(separation, bounds=(0, 5), method='bounded')
    return optimum.fun * AU_KM, optimum.x


def run(probe):
    reference, offset = EVENTS[probe]
    times, jd, vectors = read_wind(probe, offset)
    craft, comet = ephemeris(probe, 'craft'), ephemeris(probe, 'comet')
    for hours in WINDOW_HOURS:
        hits = []
        for i in np.flatnonzero(np.abs(jd-reference) < 3):
            subset = vectors[np.abs(times-times[i]) <= hours*1800]
            velocity = np.nanmedian(subset, axis=0)
            if np.all(np.isfinite(velocity)):
                dist, lag = miss(craft, comet, jd[i], velocity)
                hits.append((dist, i, lag, *velocity))
        distance, i, lag, vr, vt, vn = min(hits)
        yield dict(spacecraft=probe, median_window_hours=round(hours,6),
                   receipt_utc=datetime.fromtimestamp(times[i], timezone.utc).isoformat(),
                   miss_km=round(distance), lag_days=round(lag,6),
                   Vr_kms=round(vr,3), Vt_kms=round(vt,3), Vn_kms=round(vn,3))


if __name__ == '__main__':
    results = [row for probe in EVENTS for row in run(probe)]
    with (ROOT / 'geometry.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=results[0])
        w.writeheader()
        w.writerows(results)
