"""Shift-and-stack search in four LASCO C2 high-cadence days in June 2020."""
import argparse
import csv
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import hashlib
from itertools import product
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter

ROOT = Path(__file__).resolve().parent


def fetch(entry, cache):
    path = cache / entry['fits']
    if not path.exists():
        url = f"https://lasco-www.nrl.navy.mil/lz/level_05/{entry['day']}/c2/{path.name}"
        with urlopen(url, timeout=90) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != entry['sha256']:
            raise ValueError(path)
        path.write_bytes(data)
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != entry['sha256']:
        raise ValueError(path)
    offset = next(i for i in range(0, len(data), 80) if data[i:i+8] == b'END     ')
    start = (offset + 80 + 2879) // 2880 * 2880
    raw = np.frombuffer(data, '>i2', 1024**2, offset=start).reshape(1024, 1024)
    return raw[:290, 470:1024].astype('float32') / float(entry['exp_s'])


def find_peaks(frames, times):
    frames -= np.median(frames, axis=0)
    frames /= np.maximum(.055, 1.4826 * np.median(abs(frames), axis=0))
    np.clip(frames, -5, 12, out=frames)
    frames -= gaussian_filter(frames, (0, 4, 4))
    frames = gaussian_filter(frames, (0, 1.15, 1.15))
    for mid in np.arange(1.45, 23, 1.2):
        indices = np.flatnonzero(abs(times-mid) <= 1.40)
        if len(indices) < 10:
            continue
        proposals = []
        for vx, vy in product(range(-12, 13, 4), [32, 37, 42, 47, 52]):
            sums = np.zeros((2, 115, 500), np.float32)
            for k, i in enumerate(indices):
                dx, dy = np.rint((times[i]-mid) * np.array([vx, vy])).astype(int)
                sums[k % 2] += frames[i, 75+dy:190+dy, 35+dx:535+dx]
            n = len(indices)
            odd, even = sums / np.array([(n+1)//2, n//2])[:, None, None]
            mean, lesser = (odd+even)/2, np.minimum(odd, even)
            yy, xx = np.where((mean == maximum_filter(mean, 7)) & (lesser > .25))
            strongest = np.argsort(mean[yy, xx])[-9:]
            for y, x in zip(yy[strongest], xx[strongest]):
                proposals.append((mean[y, x], lesser[y, x], mid, x+505, y+75, vx, vy, n))
        separated = []
        for row in sorted(proposals, reverse=True):
            if all(np.hypot(row[3]-p[3], row[4]-p[4]) > 6 for p in separated):
                separated.append(row)
        yield from separated[:14]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    by_day = defaultdict(list)
    with (ROOT/'frames.csv').open() as handle:
        for entry in csv.DictReader(handle):
            by_day[entry['day']].append(entry)
    with args.output.open('w') as out:
        out.write('day,mean,smaller_half,ut_hours,column,row,vx_px_hr,vy_px_hr,frames\n')
        for day, entries in by_day.items():
            with ThreadPoolExecutor(5) as pool:
                frames = np.stack(list(pool.map(lambda e: fetch(e, args.cache), entries)))
            times = []
            for e in entries:
                h, m, s = map(int, e['utc'].split(':'))
                times.append(h + m/60 + s/3600)
            for row in find_peaks(frames, np.array(times)):
                out.write(day + ',' + ','.join(f'{v:.4f}' for v in row) + '\n')


if __name__ == '__main__':
    main()
