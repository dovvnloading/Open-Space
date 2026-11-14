"""Project the OGLE-only ephemeris onto resolved DELVE DR3 measurements."""
import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
import replay

RA, DEC = 12.57704804, -72.64897851
TARGET = '0129m727.52197'


def read_rows(name):
    with gzip.open(ROOT / 'data' / name, 'rt') as handle:
        return list(csv.DictReader(handle))


def separation(row):
    return 3600 * np.hypot((float(row['ra']) - RA) * np.cos(np.deg2rad(DEC)),
                           float(row['dec']) - DEC)


def mid_hjd(mjd, seconds):
    jd = mjd + seconds / 172800 + 2400000.5
    days = jd - 2451545
    anomaly = np.deg2rad((357.528 + .9856003 * days) % 360)
    lon = np.deg2rad((280.46 + .9856474 * days) % 360
                    + 1.915 * np.sin(anomaly) + .020 * np.sin(2*anomaly))
    obl = np.deg2rad(23.439 - .0000004 * days)
    a, b = np.deg2rad([RA, DEC])
    dot = (np.cos(lon)*np.cos(a)*np.cos(b) + np.sin(lon)*np.cos(obl)*np.sin(a)*np.cos(b)
           + np.sin(lon)*np.sin(obl)*np.sin(b))
    radius = 1.00014 - .01671*np.cos(anomaly) - .00014*np.cos(2*anomaly)
    return jd - 2450000 - radius*dot*499.005/86400


def retrieve():
    field = read_rows('field.csv.gz')
    controls = [r for r in field if separation(r) > 8 and 14 < float(r['rmag']) < 18.7]
    references = [r['objid'] for r in sorted(controls, key=separation)[:15]]
    obs = defaultdict(dict)
    for r in read_rows('measurements.csv.gz'):
        error = float(r['err'])
        if int(r['forced']) or not (0 < error < .08 and 5 < float(r['mag']) < 24):
            continue
        old = obs[r['objid']].get(r['expnum'])
        if old is None or error < float(old['err']):
            obs[r['objid']][r['expnum']] = r
    centers = {}
    for obj, exposures in obs.items():
        for band in sorted({r['filter'] for r in exposures.values()}):
            points = [[float(r[k]) for k in ('mag','ra','dec')] for r in exposures.values()
                      if r['filter'] == band]
            centers[obj, band] = np.median(points, axis=0)
    return obs, centers, references


def differential(obj, obs, centers, references):
    rows = []
    for exp, r in sorted(obs[obj].items()):
        b = r['filter']
        refs = [np.array([float(obs[k][exp][key]) for key in ('mag','ra','dec')]) - centers[k,b]
                for k in references if k != obj and exp in obs[k]]
        if len(refs) < 3:
            continue
        val = np.array([float(r[key]) for key in ('mag','ra','dec')]) - centers[obj,b]
        mag, east, north = val - np.median(refs, axis=0)
        rows.append((mid_hjd(float(r['mjd']), float(r['seconds'])), b, mag, float(r['err']),
                     east*3600*np.cos(np.deg2rad(DEC)), north*3600))
    return rows


def estimate(rows, coeff, freq, axis=2):
    if len(rows) < 12:
        return {'points': len(rows)}
    time, band, mag, err, east, north = zip(*rows)
    time = np.asarray(time)
    vals = np.array({2: mag, 4: east, 5: north}[axis])
    sigma = np.hypot(err, .010) if axis == 2 else np.full(len(rows), .035)
    template = replay.design(time, freq, 11)
    orbital = template[:,1:-2] @ coeff[1:-2]
    short = template[:,-2:] @ coeff[-2:]
    bands, ix = np.unique(band, return_inverse=True)
    dummy = np.eye(len(bands))[ix]
    X = np.column_stack((dummy, dummy*orbital[:,None], short))
    nights, group = np.unique(np.floor(time).astype(int), return_inverse=True)
    if np.linalg.matrix_rank(X) < X.shape[1] or len(nights) < 3:
        return {'points': len(rows), 'nights': len(nights)}
    w = 1/sigma
    for _ in range(5):
        c = np.linalg.lstsq(X*w[:,None], vals*w, rcond=None)[0]
        residual = vals-X@c
        w = np.sqrt(np.minimum(1, 2.5*sigma/np.maximum(abs(residual), 1e-12)))/sigma
    bread = np.linalg.pinv((X*w[:,None]).T @ (X*w[:,None]))
    by_night = np.zeros((len(nights), X.shape[1]))
    np.add.at(by_night, group, X*(w**2*residual)[:,None])
    se = np.sqrt(np.sum((by_night @ bread[:,-1])**2))
    return {'points': len(rows), 'nights': len(nights), 'scale': float(c[-1]),
            'cluster_se': float(se)}


def analyze():
    _, t, y, e, nh = next(replay.series())
    freq = json.loads((ROOT.parent/'measurements.json').read_text())['frequency_per_d']
    coeff = replay.fit(t, y, e, freq, nh)[0]
    obs, centers, refs = retrieve()
    rows = differential(TARGET, obs, centers, refs)
    subset = {b: estimate([r for r in rows if r[1] == b], coeff, freq) for b in 'griz'}
    days = sorted({int(r[0]) for r in rows})
    omitted = [estimate([r for r in rows if int(r[0]) != day], coeff, freq)['scale'] for day in days]
    field = read_rows('field.csv.gz')
    nearby = {}
    for entry in sorted(field, key=separation):
        obj = entry['objid']
        if obj == TARGET or separation(entry) > 3:
            continue
        nearby[obj] = {'arcsec': float(separation(entry)), 'r_mag': float(entry['rmag']),
                       'forced_phase': estimate(differential(obj, obs, centers, refs), coeff, freq)}
    return {'on_target': estimate(rows, coeff, freq), 'by_band': subset,
            'leave_one_night_out_range': [min(omitted),max(omitted)],
            'centroid_arcsec_per_mag': {axis: estimate(rows, coeff, freq, j)
                                        for j, axis in [(4,'east'),(5,'north')]},
            'neighbors': nearby, 'references': refs}


if __name__ == '__main__':
    print(json.dumps(analyze(), indent=2))
