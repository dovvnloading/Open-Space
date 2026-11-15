"""Replay the discovery screen from the release catalogue and per-star spectra.

Usage: python screen.py CACHE_DIRECTORY [--workers 4]
"""
import argparse
import csv
import hashlib
import itertools
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np

from analyze import CENTERS, C, ROOT, write_csv
from gzip import open as gzopen

ARCHIVE = 'https://cygnus.astro.warwick.ac.uk/phsdaj/WhiteDwarfs/'
CATALOG = 'DESI_DR1_catalogue_1.0.csv.gz'


def score(data, center, speed, width):
    w, y, _, ivar, mask = data.T
    delta = w - center * (1 + speed / C)
    good = (mask == 0) & (ivar > 0) & np.isfinite(y)
    side = good & (abs(delta) > 22) & (abs(delta) < 50)
    win = good & (abs(delta) < center * width / C)
    if side.sum() < 12 or win.sum() < 4:
        return np.nan
    design = np.column_stack((np.ones(len(w)), delta))
    xs = design[side]
    weighted = xs.T * ivar[side]
    cov = np.linalg.inv(weighted @ xs)
    coeff = cov @ weighted @ y[side]
    weights = ivar[win] / ivar[win].sum()
    projected = weights @ design[win]
    var = 1 / ivar[win].sum() + projected @ cov @ projected
    return (weights @ y[win] - projected @ coeff) / np.sqrt(var)


def inspect(options, name):
    filename = name + '.tar.gz'
    path = options.cache / filename
    url = ARCHIVE + name[2:7] + '/' + name + '/' + name + '_DESI_DR1_spectra.tar.gz'
    if not path.exists():
        tmp = path.with_suffix('.part')
        urlretrieve(url, tmp)
        tmp.replace(path)
    with tarfile.open(path) as tar:
        data = np.loadtxt(tar.extractfile(name + '_DESI_DR1_coadd-z.dat'))
    scores = []
    for width, speed in itertools.product((210, 420, 650), (-120, 0, 120)):
        z = [score(data, center, speed, width) for center in CENTERS]
        scores.append((np.sum(z) / np.sqrt(3), width, speed, z))
    available = [candidate for candidate in scores if np.isfinite(candidate[0])]
    highest = max(available, key=lambda x: x[0]) if available else (np.nan, 0, 0, [np.nan]*3)
    combined, width, speed, z = highest
    with path.open('rb') as f:
        sha = hashlib.file_digest(f, 'sha256').hexdigest()
    return dict(star=name, sha256=sha, score=combined, min_z=min(z), half_width_kms=width,
                center_shift_kms=speed, z_8500=z[0], z_8544=z[1], z_8665=z[2])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    catalog = args.cache / CATALOG
    if not catalog.exists():
        urlretrieve(ARCHIVE + 'Tables/' + CATALOG, catalog)
    with gzopen(catalog, 'rt') as f:
        rows = [r for r in csv.DictReader(f) if 'Z' in r['specType']
                and float(r['TEFF'] or 'nan') > 10000 and float(r['SN_Z'] or 'nan') > 4]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(lambda row: inspect(args, row['WDJname']), rows))
    write_csv(ROOT / 'screen.csv', sorted(results, key=lambda x: x['star']))


if __name__ == '__main__':
    main()
