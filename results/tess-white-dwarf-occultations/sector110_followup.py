"""Replay the partial S110 raw-FFI visit for the second overlapping clean S109 WD."""
import argparse
import csv
import gzip
import hashlib
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.special import ndtr

from ffi_acquire import propagated, tsv_rows
from ffi_photometry import deposit, gaia_t, match_tic
from latest_ffi import BASE, ROSTER, one
from wd_search import search

HERE = Path(__file__).resolve().parent
SOURCE = '3186021141200137472'
TIC = '55997020'
LAST = 'tess2026310211905-s0110-1-4-0312-s_ffic.fits'
FIRST = 'tess2026304234545-s0110-1-4-0312-s_ffic.fits'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def retrieve(path, url):
    if not path.exists():
        with urlopen(url, timeout=180) as response:
            path.write_bytes(response.read())


def metrics(row, coords, ratios):
    h = fits.Header.fromstring(row['wcs'], sep='\n')
    x, y = WCS(h, naxis=2).all_world2pix(*coords.T, 0)
    # First coordinate is the propagated WD, followed by its Gaia neighbours.
    if np.hypot(x[0]-row['x'], y[0]-row['y']) > 1e-6:
        raise ValueError('Cannot reproduce aperture centre')
    ox, oy = row['stamp_origin']
    yy, xx = np.indices((9, 9))
    xx, yy = xx+ox, yy+oy
    r2 = (xx-x[0])**2 + (yy-y[0])**2
    aperture = r2 <= 1.5**2
    pixels = np.array(row['stamp'], dtype=float)
    assert aperture.sum() == row['aperture_pixels']
    if not np.isclose(np.sum((pixels-row['sky'])[aperture]), row['flux'], atol=.002):
        raise ValueError('Cannot reproduce aperture photometry')
    fractions = []
    for sigma in (.65, .85, 1.05):
        # deposit sums over its entire input, so integrate one axis of the
        # neighbour-by-aperture matrix explicitly for a per-source throughput.
        dx = xx[aperture] - x[:, None]
        dy = yy[aperture] - y[:, None]
        throughput = np.sum((ndtr((dx+.5)/sigma) - ndtr((dx-.5)/sigma)) *
                            (ndtr((dy+.5)/sigma) - ndtr((dy-.5)/sigma)), axis=1)
        assert np.isclose(throughput[0], deposit(xx[aperture], yy[aperture], x[0], y[0], sigma))
        fractions.extend(1 / (1 + ratios @ throughput[1:] / throughput[0]))
    weights = np.where(r2 < 9, np.maximum(0., pixels-row['sky']), 0.)
    return dict(target_fraction=fractions[2], min_target_fraction=min(fractions),
                max_target_fraction=max(fractions), centroid_dx=float(np.sum((xx-x[0])*weights)/weights.sum()),
                centroid_dy=float(np.sum((yy-y[0])*weights)/weights.sum()))


def run(cache, out, workers):
    cache.mkdir(exist_ok=True, parents=True)
    out.mkdir(exist_ok=True, parents=True)
    stars = list(csv.DictReader((HERE/'sector110_footprint.csv').open()))
    star = next(r for r in stars if r['GaiaEDR3'] == SOURCE)
    assert (star['camera'], star['ccd']) == ('1', '4')
    queries = json.loads((HERE/'record/ffi_archive_receipts.json').read_text())
    receipts = {}
    for kind in ('gaia.tsv', 'tic82.tsv'):
        name = f'{SOURCE}-{kind}'
        path = cache/name
        url = queries[name]['url']
        retrieve(path, url)
        receipts[name] = dict(url=url, sha256=sha(path), content=path.read_text())
    tic, sep = match_tic(star, tsv_rows(cache/f'{SOURCE}-tic82.tsv'))
    assert tic and tic['TIC'] == TIC
    cone = tsv_rows(cache/f'{SOURCE}-gaia.tsv')
    assert sum(r['Source'] == SOURCE for r in cone) == 1 and all(r['Gmag'] for r in cone)
    others = [r for r in cone if r['Source'] != SOURCE]
    coords = np.array([propagated(star, 2026.8)] + [propagated(r, 2026.8) for r in others])
    ratios = np.array([[10**(-.4*(gaia_t(r, worst)-float(tic['Tmag'])-.3*worst)) for r in others]
                       for worst in (False, True)])
    roster = cache/'roster.sh'
    retrieve(roster, ROSTER)
    names = sorted(n for n in set(re.findall(r'-o (tess\d+-s0110-1-4-\d+-s_ffic\.fits)', roster.read_text()))
                   if FIRST <= n <= LAST)
    if len(names) != 2549 or names[0] != FIRST or names[-1] != LAST:
        raise ValueError('Frozen detector roster changed')
    receipts['roster'] = dict(url=ROSTER, sha256=sha(roster), first=FIRST, last=LAST, count=len(names),
                              detector_filenames_sha256=hashlib.sha256(('\n'.join(names)+'\n').encode()).hexdigest())
    receipts['fits'] = dict(base=BASE+'/s0110/{year}/{day}/1-4/{filename}', header_range=[0, 25919],
                             strip_range_fields=['start','end'], bytes_inclusive=True)
    (out/'inputs.json').write_text(json.dumps(receipts, indent=2)+'\n')
    frames = cache/'cache'
    frames.mkdir(exist_ok=True)

    def cached(n):
        path = frames/(n+'.json')
        if path.is_file():
            row = json.loads(path.read_text())
            if 'stamp' in row or row.get('status') == 'target off stampable image':
                return row
        row = one(n, star, keep_pixels=True)
        path.write_text(json.dumps(row)+'\n')
        return row

    with ThreadPoolExecutor(workers) as pool:
        rows = list(pool.map(cached, names))
    if any('error' in r for r in rows):
        raise RuntimeError('Incomplete range retrieval. Re-run with the same cache.')
    for row in rows:
        if row['flux'] is not None and row['quality'] == 0:
            row.update(metrics(row, coords, ratios))
    # Save individual input hashes and measured fields for a compact offline search replay.
    fields = [{k:v for k,v in r.items() if k not in ('stamp','stamp_origin','wcs')} for r in rows]
    payload = ''.join(json.dumps(r, separators=(',',':'), allow_nan=False)+'\n' for r in fields).encode()
    series = out/'frames.jsonl.gz'
    series.write_bytes(gzip.compress(payload,mtime=0))
    good = [r for r in fields if r.get('min_target_fraction',0) >= .80]
    t = np.array([r['time'] for r in good])
    f = np.array([r['flux'] for r in good])
    assert np.all(np.diff(t) > 0) and np.isclose(np.median(np.diff(t))*86400,200.,atol=.5)
    summary = dict(gaia=SOURCE, tic=TIC, catalogue_separation_arcsec=sep, sector=110,
                   camera=1, ccd=4, gaia_cone_rows=len(cone), attempted=len(names), headers_read=len(fields),
                   cutouts_read=sum('strip_sha256' in r for r in fields),
                   quality0_finite=sum('min_target_fraction' in r for r in fields), crowding_pass=len(good),
                   first_good_date_obs=good[0]['utc'], last_good_date_obs=good[-1]['utc'],
                   first_btjd=t[0], last_btjd=t[-1], median_flux=float(np.median(f)),
                   min_target_fraction=min(r['min_target_fraction'] for r in good),
                   nominal_target_fraction_median=float(np.median([r['target_fraction'] for r in good])),
                   frames_sha256=sha(series), search=search(t,f,200))
    (out/'search.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='search'}), flush=True)
    print('usable/periodic/isolated:', summary['search']['usable_cadences'],
          *(len(summary['search'][k]) for k in ('periodic','isolated')))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cache', type=Path)
    p.add_argument('out', type=Path)
    p.add_argument('--workers', type=int, default=9)
    a = p.parse_args()
    run(a.cache, a.out, a.workers)
