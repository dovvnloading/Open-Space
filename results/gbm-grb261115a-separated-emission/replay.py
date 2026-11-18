"""Independent frozen CTIME screen. Example: python replay.py RAW INPUTS_JSON OUTPUT_DIR"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.stats import binom

p = argparse.ArgumentParser()
p.add_argument('raw', type=Path)
p.add_argument('manifest', type=Path)
p.add_argument('output', type=Path)
p.add_argument('--compare', type=Path, help='Primary uncompressed windows.csv')
args = p.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
manifest = json.loads(args.manifest.read_text(encoding='utf-8-sig'))
for name, expected in manifest['files'].items():
    actual = hashlib.sha256((args.raw / name).read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f'hash mismatch: {name}')

records = {}
for detector in (4, 5, 6, 7, 9):
    path = args.raw / f'glg_ctime_n{detector}_bn261115924_v00.pha'
    with fits.open(path) as h:
        data = h['SPECTRUM'].data
        zero = h[0].header['TRIGTIME']
        start = np.array(data['TIME']) - zero
        end = np.array(data['ENDTIME']) - zero
        exposure = np.array(data['EXPOSURE'], dtype=float)
        counts = np.asarray(data['COUNTS'], dtype=int)[:, 2:5].sum(axis=1)
        duration = end - start
        gti = h['GTI'].data
        inside = np.logical_or.reduce([(start >= a-zero-1e-6) & (end <= b-zero+1e-6)
                                       for a, b in zip(gti['START'], gti['STOP'])])
        good = ((duration > 0) & (data['QUALITY'] == 0) & inside
                & (exposure >= .9*duration) & (exposure <= duration+.0001))
        records[detector] = (start, end, exposure, counts, good)

start, end = records[9][:2]
for detector, (s, e, _, _, _) in records.items():
    if not np.array_equal(s, start) or not np.array_equal(e, end):
        raise ValueError(f'detector {detector} has a different time grid')

sidebands = {}
for detector, (s, e, exposure, counts, good) in records.items():
    sidebands[detector] = []
    for low, high in [(-72, -36), (44, 72)]:
        take = good & (s >= low) & (e <= high)
        live = float(exposure[take].sum())
        if live < .9*(high-low):
            raise ValueError(f'insufficient sideband exposure {detector}')
        sidebands[detector].append((int(counts[take].sum()), live))

def test(detector, indices):
    _, _, exp, counts, _ = records[detector]
    observed, live = int(counts[indices].sum()), float(exp[indices].sum())
    p = max(float(binom.sf(observed-1, observed+b, live/(live+bg_live)))
            for b, bg_live in sidebands[detector])
    return p, observed, live

regions = [('pre',-32,-3), ('late',8,40), ('early_control',-106,-77), ('late_control',77,109)]
all_results = {}
for name, low, high in regions:
    windows = []
    for width in (.256, .512, 1.024, 2.048):
        for i in np.flatnonzero((start >= low) & (end <= high)):
            j = int(np.searchsorted(end, start[i] + width - 1e-5))
            if j >= len(start) or end[j] > high or abs(end[j]-start[i]-width) > 1e-4:
                continue
            ids = np.arange(i,j+1)
            dt = end[ids]-start[ids]
            if (i == 0 or j+1 == len(start) or np.max(abs(dt-dt[0])) > 1e-4
                or np.max(abs(start[ids[1:]]-end[ids[:-1]]), initial=0) > 1e-4
                or abs(end[i-1]-start[i-1]-dt[0]) > 1e-4
                or abs(end[j+1]-start[j+1]-dt[0]) > 1e-4
                or abs(start[i]-end[i-1]) > 1e-4 or abs(start[j+1]-end[j]) > 1e-4
                or not all(np.all(r[4][ids]) for r in records.values())):
                continue
            checks = {str(d): test(d, ids) for d in (9,6,7,4,5)}
            windows.append({'start':float(start[i]), 'end':float(end[j]), 'width':width,
                            'joint_p':max(checks['9'][0], checks['6'][0]), 'detectors':checks})
    windows.sort(key=lambda x:x['joint_p'])
    all_results[name] = windows
ntrials = len(all_results['pre']) + len(all_results['late'])
prompt = {}
for detector, (s,e,ex,c,good) in records.items():
    ids = np.flatnonzero(good & (s >= -1) & (e <= 4))
    prompt[str(detector)] = {'bins':len(ids),'p_counts_live':test(detector, ids)}
summary = {'trials':ntrials, 'regions':{}, 'prompt':prompt, 'sidebands':sidebands}
for name, _, _ in regions:
    w = all_results[name]
    summary['regions'][name] = {'count':len(w), 'best':w[0] if w else None,
                                'crossings':sum(x['joint_p']*ntrials <= .01 for x in w)}
if args.compare:
    peer_windows = {
        (region.replace('early', 'pre'), row['start'], row['end'], row['width']): row
        for region, rows in all_results.items()
        for row in rows
    }
    comparisons = 0
    with args.compare.open(encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream):
            key = (row['region'], float(row['start_s']), float(row['stop_s']),
                   float(row['nominal_width_s']))
            window = peer_windows.pop(key)
            if window['joint_p'] != float(row['joint_p']):
                raise ValueError(f'Joint score mismatch for {key}')
            for detector, triplet in window['detectors'].items():
                primary = (float(row[f'n{detector}_p']), int(row[f'n{detector}_counts']),
                           float(row[f'n{detector}_exposure_s']))
                if tuple(triplet) != primary:
                    raise ValueError(f'Detector {detector} mismatch for {key}')
            comparisons += 1
    if peer_windows:
        raise ValueError(f'{len(peer_windows)} unmatched independent windows')
    summary['compared_windows'] = comparisons
with (args.output / 'replay.json').open('x', encoding='utf-8', newline='\n') as stream:
    stream.write(json.dumps(summary, indent=2, allow_nan=False) + '\n')
print(json.dumps(summary, indent=2))
