"""Time held-out TESS primary eclipses. Fetch the immutable manifest into a local cache."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.signal import argrelmin
from scipy.stats import binned_statistic

from reproduce import ROOT, EPOCH, closest, orbit, separation


def lightcurve(data):
    """Read the first noncompressed scalar SPOC binary table; reject other layouts."""
    start = data.index(b'XTENSION')
    header = {}
    for off in range(start, len(data), 80):
        card = data[off:off+80].decode('ascii')
        if card.startswith('END '):
            break
        if card[8:10] == '= ':
            header[card[:8].strip()] = card[10:].split('/')[0].strip().strip("'").strip()
    fmt = {'D': '>f8', 'E': '>f4', 'J': '>i4'}
    cols = [(header[f'TTYPE{i}'], fmt[header[f'TFORM{i}']])
            for i in range(1, int(header['TFIELDS'])+1)]
    dtype = np.dtype(cols)
    assert header['XTENSION'] == 'BINTABLE' and dtype.itemsize == int(header['NAXIS1'])
    assert float(header['BJDREFI']) + float(header['BJDREFF']) == 2457000
    return np.frombuffer(data, dtype=dtype, count=int(header['NAXIS2']), offset=(off//2880+1)*2880)


def eclipse(time, flux, predicted):
    near = np.abs(time - predicted) < .34
    x, y = time[near]-predicted, flux[near].copy()
    out = np.abs(x) > .14
    if len(x) < 200 or out.sum() < 45 or (np.abs(x)<.05).sum() < 35:
        return None
    y /= np.polyval(np.polyfit(x[out], y[out], 1), x)
    edges = np.arange(-.12, .121, .005)
    b = binned_statistic(x, y, statistic='median', bins=edges).statistic
    centers = (edges[1:]+edges[:-1])/2
    depth = 1-np.nanmedian(b[np.abs(centers)<.018])
    if not np.isfinite(depth) or depth < .02:
        return None
    level = 1-depth/2
    changes = np.where(np.diff(b < level))[0]
    changes = [j for j in changes if np.isfinite(b[j:j+2]).all()]
    left = [j for j in changes if -.1 < centers[j] < 0]
    right = [j for j in changes if 0 < centers[j] < .1]
    if len(left) != 1 or len(right) != 1:
        return None
    def crossing(j):
        return centers[j] + (level-b[j])*(centers[j+1]-centers[j])/(b[j+1]-b[j])
    return (crossing(left[0])+crossing(right[0]))/2, depth


def validate(cache):
    apparent, radii = orbit(5700)
    grid = np.arange(3500, 5650, .12)
    pos = apparent(grid)
    d = separation(pos, 0, 1)
    ix = argrelmin(d)[0]
    ix = ix[(d[ix]<radii[0]) & (pos[ix, 1, 2]>pos[ix, 0, 2])]
    primary = [closest(apparent, 0, 1, grid[k]-.12, grid[k]+.12) for k in ix]
    observations, sectors = [], []
    with (ROOT/'tess-products.csv').open() as fh:
        products = list(csv.DictReader(fh))
    for row in sorted(products, key=lambda r: int(r['sector'])):
        sec = int(row['sector'])
        filename = cache / f'tess-{sec}.fits'
        if not filename.exists():
            url = 'https://mast.stsci.edu/api/v0.1/Download/file?' + urlencode({'uri': row['mast_uri']})
            filename.write_bytes(urlopen(url, timeout=120).read())
        data = filename.read_bytes()
        if hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ValueError(f'Sector {sec} differs from frozen product')
        table = lightcurve(data)
        good = (table['QUALITY']==0) & np.isfinite(table['TIME']) & np.isfinite(table['PDCSAP_FLUX'])
        time = table['TIME'][good] + (2457000 - EPOCH)
        flux = table['PDCSAP_FLUX'][good]
        typical_error = float(np.median(table['PDCSAP_FLUX_ERR'][good])/np.median(flux))
        flux = flux/np.median(flux)
        eligible = [t for t in primary if time.min()+.4 < t < time.max()-.4]
        sectors.append({'sector': sec, 'eligible_eclipses': len(eligible), 'good_cadences': len(time),
                        'start_bjd': round(EPOCH+min(time), 5), 'end_bjd': round(EPOCH+max(time), 5),
                        'median_fractional_cadence_error': round(typical_error, 5)})
        for t in eligible:
            result = eclipse(time, flux, t)
            if result:
                offset, depth = result
                observations.append({'sector':sec, 'model_bjd': round(EPOCH+t, 6),
                                     'observed_minus_model_minutes': round(offset*1440, 2),
                                     'fractional_depth': round(depth, 4)})
    return {'sectors':sectors, 'observations':observations}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    print(json.dumps(validate(args.cache), indent=2))
