"""Replay the broad, two-body JUICE preselection with the JPL element tables."""
import argparse
import csv
import gzip
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline

EPOCH = 2461200.5
LIMIT = 0.012


def read_orbits(directory):
    names, elements = [], []
    for label, start in [('NUMBR', 25), ('UNNUM', 14)]:
        with gzip.open(directory / f'ELEMENTS.{label}.gz', 'rt') as rows:
            next(rows)
            next(rows)
            for row in rows:
                fields = row[start:].split()
                if fields[0] == '61200' and 1.3 < float(fields[1]) < 5.5:
                    names.append(row[:start].strip())
                    elements.append([float(value) for value in fields[1:8]])
    return np.array(names), np.array(elements).T


def candidates(names, elements, craft):
    a, ecc, inc, peri, node, anomaly, magnitude = elements
    inc, peri, node, anomaly = np.radians([inc, peri, node, anomaly])
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(peri), np.sin(peri)
    co, so = np.cos(node), np.sin(node)
    p = np.column_stack([co*cw-so*sw*ci, so*cw+co*sw*ci, sw*si])
    q = np.column_stack([-co*sw-so*cw*ci, -so*sw+co*cw*ci, cw*si])
    normal = np.cross(p, q)
    mean_motion = 0.01720209895 / a**1.5
    daily = CubicHermiteSpline(craft[:, 0], craft[:, 1:4], craft[:, 4:], axis=0)
    radii_min = a * (1-ecc)
    radii_max = a * (1+ecc)
    best = np.full(len(a), LIMIT)
    dates = np.full(len(a), np.nan)
    for start in np.arange(craft[0, 0], craft[-1, 0], 4.0):
        half = min(4.0, craft[-1, 0] - start) / 2
        jd = start + half
        pos, vel = daily(jd), daily(jd, 1)
        radius = np.linalg.norm(pos)
        # The radial and plane gates pad the interval for spacecraft motion.
        mask = (radius >= radii_min-.08) & (radius <= radii_max+.08)
        mask &= abs(normal @ pos) <= best+.065
        idx = np.flatnonzero(mask)
        e = ecc[idx]
        mean = anomaly[idx] + (jd-EPOCH)*mean_motion[idx]
        E = mean + e*np.sin(mean)
        for _ in range(8):
            E -= (E-e*np.sin(E)-mean)/(1-e*np.cos(E))
        c, s = np.cos(E), np.sin(E)
        beta = np.sqrt(1-e**2)
        plane_x, plane_y = p[idx], q[idx]
        r = a[idx,None] * ((c-e)[:,None]*plane_x + (beta*s)[:,None]*plane_y) - pos
        v = (a[idx]*mean_motion[idx]/(1-e*c))[:,None] * (
            -s[:,None]*plane_x + (beta*c)[:,None]*plane_y) - vel
        dt = np.clip(-np.sum(r*v, axis=1)/np.sum(v*v, axis=1), -half, half)
        miss = np.linalg.norm(r+dt[:,None]*v, axis=1)
        keep = miss < best[idx]
        best[idx[keep]] = miss[keep]
        dates[idx[keep]] = jd + dt[keep]
    for j in np.argsort(best):
        if best[j] >= LIMIT:
            break
        yield names[j], magnitude[j], round(best[j], 9), round(dates[j], 7)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('catalog_directory', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    daily = np.loadtxt(Path(__file__).parent/'data/juice-daily.csv.gz', delimiter=',', skiprows=1)
    names, elements = read_orbits(args.catalog_directory)
    with args.output.open('w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['name','H_mag','screen_au','screen_jd_tdb'])
        writer.writerows(candidates(names, elements, daily))
    print(f'{len(names):,} orbits searched')
