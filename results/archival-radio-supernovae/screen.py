"""Repeat the initial 3 GHz fixed-position measurements for the frozen BTS sample."""
import argparse
import hashlib
import io
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parent
CUTOUT = 'https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync?'


def fit(data, row):
    with fits.open(io.BytesIO(data)) as hdus:
        header = hdus[0].header
        wcs = WCS(header).celestial
        image = 1000 * hdus[0].data.squeeze().astype(float)
    if not all(np.isclose(abs(header[key]) * 3600, 1, atol=.01)
               for key in ('CDELT1', 'CDELT2')):
        raise ValueError('Expected one-arcsecond pixels')
    x, y = map(float, wcs.world_to_pixel_values(row.ra, row.dec))
    yy, xx = np.indices(image.shape)
    distance = np.hypot(xx-x, yy-y)
    background = image[(distance >= 8) & (distance < 30) & np.isfinite(image)]
    core = image[distance < 2]
    if len(background) < 650 or len(core) < 8 or not np.isfinite(core).all():
        raise ValueError('Incomplete cutout')
    median = np.median(background)
    rms = 1.4826 * np.median(abs(background - median))
    if not rms > 0:
        raise ValueError('Invalid field noise')
    image -= median
    iy, ix = np.unravel_index(np.argmax(np.where(distance <= 2, image, -np.inf)), image.shape)
    return dict(
        date=header['DATE-OBS'], mjd=float(Time(header['DATE-OBS']).mjd),
        flux=float(map_coordinates(image, [[y], [x]], order=1)[0]),
        peak=image[iy, ix], rms=rms, offset=distance[iy, ix],
        bmaj=header['BMAJ'] * 3600, bmin=header['BMIN'] * 3600,
        sha256=hashlib.sha256(data).hexdigest(),
    )


def cutout_url(row, radius):
    circle = f'{row.ra:.7f} {row.dec:.7f} {radius:.7f}'
    return CUTOUT + urlencode({'ID': row.uri, 'CIRCLE': circle})


def download(row, radius, path):
    if not path.exists():
        with urlopen(cutout_url(row, radius), timeout=90) as response:
            data = response.read()
        if not data.startswith(b'SIMPLE'):
            raise ValueError('Archive did not return FITS')
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    return path.read_bytes()


def repeat(row, cache):
    result = dict(name=row.name, ztf=row.ztf, epoch=row.epoch)
    try:
        data = download(row, .009, cache / f'{row.ztf}_{row.epoch}.fits')
        result.update(fit(data, row))
    except (OSError, ValueError) as exc:
        result['error'] = str(exc)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    rows = pd.read_csv(ROOT / 'screen.csv.gz')
    args.cache.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.workers) as workers:
        output = workers.map(lambda row: repeat(row, args.cache), rows.itertuples(index=False))
        pd.DataFrame(output).to_csv(args.cache / 'replayed.csv', index=False)
