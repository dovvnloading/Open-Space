"""Refetch and hash original products. Requires NumPy and Astropy for FITS."""
import argparse
import csv
import gzip
import hashlib
import json
import shutil
import tarfile
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits

from analyze import BAND, ROOT, fit, write_csv

SDSS = 'https://data.sdss.org/sas/dr17/'
WARWICK = 'https://cygnus.astro.warwick.ac.uk/phsdaj/WhiteDwarfs/'
STAR = 'WDJ230752.92-000207.45'
SOURCES = {
    'SDSS-2000': ('spec-0381-51811-0425.fits', SDSS + 'sdss/spectro/redux/26/spectra/lite/0381/spec-0381-51811-0425.fits'),
    'BOSS-2010': ('spec-4209-55478-0518.fits', SDSS + 'eboss/spectro/redux/v5_13_2/spectra/lite/4209/spec-4209-55478-0518.fits'),
    'LAMOST-2015': ('lamost-dr5.fits.gz', 'https://www.lamost.org/dr5/v3/spectrum/fits/367608207'),
    'DESI-2021': (STAR + '.tar.gz', WARWICK + 'J2307/' + STAR + '/' + STAR + '_DESI_DR1_spectra.tar.gz'),
    'SDSS-plate': ('plate0381.fits', SDSS + 'sdss/spectro/redux/26/0381/spPlate-0381-51811.fits'),
    'BOSS-plate': ('plate4209.fits', SDSS + 'eboss/spectro/redux/v5_13_2/4209/spPlate-4209-55478.fits')
}


def read_spectrum(epoch, path):
    if epoch == 'DESI-2021':
        with tarfile.open(path) as tar:
            cols = np.loadtxt(tar.extractfile(STAR + '_DESI_DR1_coadd-z.dat')).T
        return cols[0], cols[1], cols[3], cols[4], cols[4]
    if epoch == 'LAMOST-2015':
        flux, ivar, wave, andmask, ormask = fits.getdata(path)
    else:
        t = fits.getdata(path, 1)
        wave = 10 ** t['loglam'].astype(float)
        flux, ivar, andmask, ormask = (t[name] for name in ('flux', 'ivar', 'and_mask', 'or_mask'))
    return wave, flux, ivar, andmask, ormask


def controls(path, epoch):
    rows = []
    with fits.open(path) as h:
        flux, ivar, mask = (h[i].data for i in range(3))
        hdr = h[0].header
        wave = 10**(hdr['COEFF0'] + np.arange(flux.shape[1]) * hdr['COEFF1'])
        for i, plug in enumerate(h[5].data):
            res = fit(wave, flux[i], ivar[i], mask[i])
            if res is not None:
                rows.append(dict(epoch=epoch, fiber=int(plug['FIBERID']),
                                 plug_type=plug['OBJTYPE'].strip(), **res))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path, help='directory for six survey products (~165 MiB)')
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    receipts = []
    control_rows = []
    spec_rows = []
    for epoch, (filename, url) in SOURCES.items():
        p = args.cache / filename
        if not p.exists():
            tmp = p.with_suffix('.part')
            with urlopen(url, timeout=120) as response, tmp.open('wb') as fp:
                shutil.copyfileobj(response, fp)
            tmp.replace(p)
        with p.open('rb') as handle:
            digest = hashlib.file_digest(handle, 'sha256').hexdigest()
        receipts.append(dict(epoch=epoch, file=filename, url=url, bytes=p.stat().st_size, sha256=digest))
        if epoch.endswith('plate'):
            control_rows.extend(controls(p, epoch))
        else:
            w, y, iv, am, om = read_spectrum(epoch, p)
            selected = (w > 8300) & (w < 8850)
            for values in zip(w[selected], y[selected], iv[selected], am[selected], om[selected]):
                w0, y0, iv0, a0, o0 = values
                spec_rows.append((epoch, float(w0), float(y0), float(iv0), int(a0), int(o0)))
    with gzip.open(ROOT / 'spectra.csv.gz', 'wt', newline='') as fp:
        writer = csv.writer(fp)
        writer.writerow(['epoch', 'wave', 'flux', 'ivar', 'andmask', 'ormask'])
        writer.writerows(spec_rows)
    with gzip.open(ROOT / 'plate_controls.csv.gz', 'wt', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=control_rows[0])
        writer.writeheader()
        writer.writerows(control_rows)
    write_csv(ROOT / 'sources.csv', receipts)


if __name__ == '__main__':
    main()
