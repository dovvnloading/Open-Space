"""Download small tangent-plane cutouts of GO-14618 calibrated MAST products.

The CTE-corrected FLC frames are in electrons, DRC frames in electrons/s. For
FLC use the chip header's SIP WCS without its small lookup-table corrections.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import io
from pathlib import Path
import time
from urllib.request import Request, urlopen

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from scipy.ndimage import map_coordinates

RA = 187.68775
DEC = 12 + 23 / 60 + 26.44 / 3600
URL = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/'
SPLITS = ['id5o42lhq', 'id5o42ljq', 'id5o43m6q', 'id5o43m8q', 'id5o52m9q', 'id5o52mbq',
          'id5o43m1q', 'id5o43m2q', 'id5o43m4q']


def fetch_range(filename, start, stop):
    request = Request(URL + filename, headers={'Range': f'bytes={start}-{stop - 1}'})
    for attempt in range(3):
        try:
            with urlopen(request, timeout=60) as response:
                content_range = response.headers.get('Content-Range', '')
                if not content_range.startswith(f'bytes {start}-{stop - 1}/'):
                    raise ValueError(f'Unexpected byte range {content_range}')
                return response.read()
        except (OSError, ValueError):
            if attempt == 2:
                raise
            time.sleep(attempt + 1)


def download(filename):
    split = filename.endswith('_flc.fits')
    output = Path('stamps') / (filename.replace('.fits', '') + '.npz' if split else filename + '.npz')
    if output.exists():
        return
    offset = 0
    for hdu_index in range(8):
        size = 57600
        while True:
            block = io.BytesIO(fetch_range(filename, offset, offset + size))
            try:
                header = fits.Header.fromfile(block)
                break
            except OSError:
                if size >= 460800:
                    raise
                size *= 2
        start = offset + block.tell()
        shape = [header[f'NAXIS{i+1}'] for i in range(header['NAXIS'])]
        size = int(np.prod(shape)) * (abs(header['BITPIX']) // 8) if shape else 0
        offset = start + ((size + 2879) // 2880) * 2880
        if hdu_index == 0:
            primary = header
            continue
        if header.get('EXTNAME') != 'SCI':
            continue
        if split:
            for key in list(header):
                if key.startswith(('CPDIS', 'DP1', 'DP2', 'D2IM')):
                    del header[key]
        wcs = WCS(header, naxis=2)
        nx, ny = shape
        x, y = wcs.world_to_pixel_values(RA, DEC)
        if not (50 < x < nx - 50 and 50 < y < ny - 50):
            continue
        if header['BITPIX'] != -32:
            raise ValueError('Expected float32 SCI array')
        yy, xx = np.mgrid[-48:49, -48:49]
        sky_ra = RA + .04 * xx / 3600 / np.cos(np.deg2rad(DEC))
        gx, gy = wcs.world_to_pixel_values(sky_ra, DEC + .04 * yy / 3600)
        lo = int(np.floor(gy.min())) - 2
        hi = int(np.ceil(gy.max())) + 3
        pixels = fetch_range(filename, start + lo * nx * 4, start + hi * nx * 4)
        strip = np.frombuffer(pixels, dtype='>f4').reshape(-1, nx)
        image = map_coordinates(strip, [gy - lo, gx], order=1, mode='constant', cval=np.nan)
        np.savez(output, image=image, filter=primary['FILTER'], mjd=primary['EXPSTART'],
                 exptime=primary['EXPTIME'], unit=header['BUNIT'])
        return
    raise ValueError(f'Source outside supported science extensions in {filename}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    Path('stamps').mkdir(exist_ok=True)
    visits = [f'id5o{visit:02d}0{sequence}0_drc.fits' for visit in range(1, 54) for sequence in (1, 2)]
    products = visits + [name + '_flc.fits' for name in SPLITS]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(download, products))
