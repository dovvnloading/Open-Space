"""Rebuild the frozen, unscaled SCI cutout payload from the individual MAST products."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
import json
from urllib.request import Request, urlopen

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np

from analyze import POSITIONS, ROOT


def read_range(name, start, end):
    url = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/' + name
    with urlopen(Request(url, headers={'Range': f'bytes={start}-{end}'}), timeout=120) as response:
        content_range = response.headers.get('Content-Range', '')
        data = response.read()
    if not content_range.startswith(f'bytes {start}-{end}/') or len(data) != end-start+1:
        raise ValueError(f'Incomplete HTTP range for {name}')
    return data


def extract(item):
    key, info = item
    number, site, *_ = key.split('-')
    name = info['product']
    buf = io.BytesIO(read_range(name, 0, 288000-1))
    fits.Header.fromfile(buf)
    h = fits.Header.fromfile(buf)
    offset = buf.tell()
    if h.get('EXTNAME') != 'SCI' or h['BITPIX'] != -32:
        raise ValueError(f'Unexpected SCI encoding in {name}')
    wcs = WCS(h, naxis=2)
    c = SkyCoord(POSITIONS[int(number)][int(site)])
    x, y = np.rint(wcs.world_to_pixel(c)).astype(int)
    w, hgt, r = h['NAXIS1'], h['NAXIS2'], 38
    if min(x-r, y-r, w-1-x-r, hgt-1-y-r) < 0:
        raise ValueError(f'Outside image for {key}')
    stripe = read_range(name, offset+4*w*(y-r), offset+4*w*(y+r+1)-1)
    a = np.frombuffer(stripe, dtype='>f4').reshape(-1, w)[:, x-r:x+r+1].copy()
    if hashlib.sha256(a.tobytes()).hexdigest() != info['sha256_pixels_be_float32']:
        raise ValueError(f'MAST SCI pixels differ from the frozen release for {key}')
    stored = WCS(fits.Header.fromstring(info['cutout_header'], sep='\n'), naxis=2)
    test = np.array([[0., 0.], [38., 38.], [76., 76.]])
    difference = stored.all_pix2world(test, 0) - wcs.all_pix2world(test + [x-r, y-r], 0)
    if np.max(np.abs(difference)) > 1e-7:
        raise ValueError(f'MAST astrometry differs for {key}')
    return key, a


if __name__ == '__main__':
    items = json.loads((ROOT / 'manifest.json').read_text()).items()
    with ThreadPoolExecutor(4) as pool:
        data = dict(pool.map(extract, items))
    np.savez_compressed(ROOT / 'cutouts.npz', **data)
    print(f'Verified {len(data)} cutouts from MAST.')
