"""Fit three quasar line profiles from the archived pixels or re-extract them from DESI."""
import argparse
import csv
import gzip
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from measure import fit_line

TARGETS = [39627745021400761, 39627745021400804, 39627745021400959, 39627745021400845]
LINES = [('CIV', 1549.48, 'B'), ('CIII', 1908.73, 'B'), ('MgII', 2798.75, 'R')]
REFERENCE = 2.511


def extract(source):
    from astropy.io import fits
    with fits.open(source) as hdus, gzip.open(HERE / 'spectra.csv.gz', 'wt') as out:
        writer = csv.writer(out)
        writer.writerow(['member', 'line', 'wave', 'flux', 'ivar'])
        for target, name in zip(TARGETS, 'ABCD'):
            row, = np.flatnonzero(hdus['FIBERMAP'].data['TARGETID'] == target)
            for line, rest, arm in LINES:
                wave = hdus[arm+'_WAVELENGTH'].data
                flux = hdus[arm+'_FLUX'].data[row]
                ivar = hdus[arm+'_IVAR'].data[row]
                mask = hdus[arm+'_MASK'].data[row]
                good = (abs(wave / REFERENCE - rest) < 110) & (ivar > 0) & (mask == 0)
                writer.writerows((name, line, w, f, iv) for w, f, iv in zip(wave[good], flux[good], ivar[good]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--extract', metavar='COADD')
    parser.add_argument('--bootstrap', type=int, default=0,
                        help='noise-only parametric bootstrap count for Mg II width intervals')
    args = parser.parse_args()
    if args.extract:
        extract(args.extract)
    rng = np.random.default_rng(101339)
    with gzip.open(HERE / 'spectra.csv.gz', 'rt') as stream:
        pixels = list(csv.DictReader(stream))
    for name in 'ABCD':
        for line, rest, arm in LINES:
            cols = [[float(p[k]) for p in pixels if p['member'] == name and p['line'] == line]
                    for k in ['wave', 'flux', 'ivar']]
            w, f, iv = np.array(cols)
            expected = rest * REFERENCE
            p, model = fit_line(w, f, iv, expected)
            width = lambda par: 2.35482 * 299792.458 * par[4] / par[3]
            print(name, line, f'z={p[3]/rest-1:.5f}, FWHM={width(p):.0f} km/s', end='')
            if line == 'MgII' and args.bootstrap:
                widths = []
                centers = []
                for _ in range(args.bootstrap):
                    flux = model(p, w) + rng.normal(size=len(w)) / np.sqrt(iv)
                    fp, _ = fit_line(w, flux, iv, expected)
                    widths.append(width(fp))
                    centers.append(fp[3] / rest - 1)
                low, high = np.percentile(widths, [2.5, 97.5])
                lo_z, hi_z = np.percentile(centers, [2.5, 97.5])
                print(f', conditional noise interval [{low:.0f}, {high:.0f}], z [{lo_z:.5f}, {hi_z:.5f}]', end='')
            print(flush=True)


if __name__ == '__main__':
    main()
