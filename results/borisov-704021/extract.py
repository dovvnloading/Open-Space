"""Extract one calibrated MegaCam cutout using the SSOS footprint position.

Usage: python extract.py /download/path/2453634p.fits.fz
The separate topocentric Horizons position in imaging-ephemeris.csv sets the photometry aperture.
"""
import argparse
import csv
from pathlib import Path
import warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

ROOT = Path(__file__).resolve().parent


def extract(path):
    image_id = path.name.split('.')[0]
    with open(ROOT / 'data/archive-images.tsv') as f:
        record = next(r for r in csv.DictReader(f, delimiter='\t') if r['Image'] == image_id)
    coords = [[float(record['Object_RA']), float(record['Object_Dec'])]]
    with fits.open(path) as hdus:
        for chip, hdu in enumerate(hdus[1:], 1):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FITSFixedWarning)
                x, y = WCS(hdu.header).all_world2pix(coords, 0)[0]
            height, width = hdu.shape
            if not (80 < x < width - 80 and 80 < y < height - 80):
                continue
            x0, y0 = int(x) - 70, int(y) - 70
            data = hdu.section[y0:y0+141, x0:x0+141].astype(np.float32)
            header = hdu.header.copy()
            header['CRPIX1'] -= x0
            header['CRPIX2'] -= y0
            header['EXTNO'] = chip
            header['RA_PRED'], header['DEC_PRED'] = coords[0]
            header['X_PRED'], header['Y_PRED'] = x-x0, y-y0
            for key in ['MJD-OBS', 'EXPTIME', 'DATE-OBS', 'FILTER']:
                header[key] = hdus[0].header.get(key, hdu.header.get(key, ''))
            fits.PrimaryHDU(data, header).writeto(ROOT / 'data' / (image_id + '.fits'), overwrite=True)
            return
    raise ValueError(f'{image_id} has no complete image stamp at the footprint position')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('image', type=Path)
    extract(parser.parse_args().image)
