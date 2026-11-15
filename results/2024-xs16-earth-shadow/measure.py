"""Compare exposure-averaged geometric Earth shadow with MPC L34 photometry."""
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent
C = 299792.458
AU = 149597870.7
RS = 695700.
DAY = 86400.
# In December 2024 TT-UTC = 69.184 s; TDB-TT is below 0.002 s.
TT_UTC = 69.184


def jd(instant):
    return instant.timestamp() / DAY + 2440587.5


def iso(tdb):
    seconds = (tdb - 2440587.5) * DAY - TT_UTC
    return datetime.fromtimestamp(seconds, timezone.utc).isoformat(timespec='milliseconds')


def vectors():
    table = np.genfromtxt(ROOT / 'vectors.csv', delimiter=',', names=True)
    cols = table.dtype.names
    def spline(i):
        pos = np.column_stack([table[k] for k in cols[i:i+3]])
        vel = np.column_stack([table[k] for k in cols[i+3:i+6]])
        return CubicHermiteSpline(table['jd_tdb'], pos, vel * DAY, extrapolate=False)
    return spline(1), spline(7)


def angles(body, sun, t, re):
    x = body(t)
    v = sun(t) - x
    a, b = np.linalg.norm(x, axis=-1), np.linalg.norm(v, axis=-1)
    cos_term = np.sum(-x*v, axis=-1)
    gap = np.arctan2(np.linalg.norm(np.cross(-x, v), axis=-1), cos_term)
    return gap, np.arcsin(re / a), np.arcsin(RS / b)


def covered(sep, earth, sun):
    if sep >= earth + sun:
        return 0.
    if sep <= abs(earth - sun):
        return min(1., (earth / sun)**2)
    # Area of intersection of two apparent disks, normalized to solar area.
    x = (sep**2 + sun**2 - earth**2) / (2 * sep * sun)
    y = (sep**2 + earth**2 - sun**2) / (2 * sep * earth)
    radicand = (-sep + sun + earth) * (sep + sun - earth) * (sep - sun + earth) * (sep + sun + earth)
    area = sun**2 * np.arccos(np.clip(x, -1, 1)) + earth**2 * np.arccos(np.clip(y, -1, 1))
    return (area - .5 * np.sqrt(max(0., radicand))) / (np.pi * sun**2)


def analysis():
    body, sun = vectors()
    with (ROOT / 'mpc.csv').open() as f:
        obs = [r for r in csv.DictReader(f) if r['stn'] == 'L34' and r['mag']]
    with (ROOT / 'l34-unshadowed.csv').open() as f:
        model = list(csv.DictReader(f))
    obsjd = [jd(datetime.strptime(r['utc'], '%Y-%b-%d %H:%M').replace(tzinfo=timezone.utc)) for r in model]
    ranges = [float(r['range_au']) for r in model]
    mags = [float(r['horizons_V']) for r in model]

    def emission(t):
        return t + TT_UTC/DAY - np.interp(t, obsjd, ranges) * AU/C/DAY

    def contact(t, interior):
        sep, re, rs = angles(body, sun, emission(t), 6371.)
        return sep - re + (rs if interior else -rs)

    grid = np.linspace(obsjd[0], obsjd[-1], 1150)
    contacts = []
    for interior in [False, True]:
        values = [contact(t, interior) for t in grid]
        for i in np.flatnonzero(np.diff(np.sign(values))):
            t = brentq(contact, grid[i], grid[i+1], args=(interior,), xtol=5e-10, rtol=1e-15)
            contacts.append(iso(t + TT_UTC/DAY))
    results = []
    for r in obs:
        t = jd(datetime.fromisoformat(r['obstime'].replace('Z', '+00:00')))
        times = t + np.linspace(-.5, .5, 221) * float(r['exp']) / DAY
        def loss(radius, offset=0.):
            params = angles(body, sun, emission(times + offset / DAY), radius)
            samples = np.array([covered(*a) for a in zip(*params)])
            return float(np.trapezoid(samples, dx=1/(len(samples)-1)))
        nominal = loss(6371.)
        l, u = loss(6356.752), loss(6378.137)
        baseline = np.interp(t, obsjd, mags)
        results.append(dict(utc=r['obstime'], observed_G=float(r['mag']), rms=float(r['rmsmag']),
                            exposure_s=int(r['exp']), mean_covered=nominal, no_shadow_V=float(baseline),
                            radius_stress_fraction=[l,u], time_stress_fraction=[loss(6371.,-30), loss(6371.,30)]))
    a, b = results[2:4]
    observed = b['observed_G'] - a['observed_G']
    baseline = b['no_shadow_V'] - a['no_shadow_V']
    darkening = -2.5*np.log10((1-b['mean_covered']) / (1-a['mean_covered']))
    return dict(contacts_utc=sorted(contacts),observations=results,observed_change_mag=observed,
                quoted_difference_error_mag=np.hypot(a['rms'], b['rms']),smooth_change_mag=baseline,
                eclipse_change_mag=darkening,combined_change_mag=baseline+darkening)


if __name__ == '__main__':
    print(json.dumps(analysis(), indent=2))
