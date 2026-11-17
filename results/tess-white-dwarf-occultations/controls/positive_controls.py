"""Blindly search every inventoried SPOC light curve for three known systems.

Use control_inventory.json from the same directory, an input cache, and the
survey's wd_search.py on PYTHONPATH. Reuses the repository's eight-column
digest defined in results/tess-additional-transits/fetch.py when in a checkout.
"""
import argparse
import concurrent.futures
import hashlib
import runpy
import json
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from wd_search import search

HERE = Path(__file__).resolve().parent
BASE = 'https://stpubdata.s3.amazonaws.com/'
# TIC: published period, epoch (BTJD), reference.
KNOWN = {
    267574918: (1.407939217, 2038.4358981, 'https://doi.org/10.1093/mnras/stad766'),
    60040774: (.404807, 2201.762996, 'https://doi.org/10.1093/mnras/stac2197'),
    219244444: (.3037036366, -5476.951433, 'https://academic.oup.com/mnras/article/376/2/919/1077105'),
}


def obtain(key, expected_size, cache):
    path = cache / key.split('/')[-1]
    if not path.exists() or path.stat().st_size != expected_size:
        for attempt in range(3):
            try:
                with urllib.request.urlopen(BASE + key, timeout=90) as response:
                    data = response.read()
                if len(data) != expected_size:
                    raise OSError(f'expected {expected_size} bytes, got {len(data)}')
                tmp = path.with_suffix('.pending')
                tmp.write_bytes(data)
                tmp.replace(path)
                break
            except OSError:
                if attempt == 2:
                    raise
    return path


def measure(row, product, cache, digest):
    tic = row['tic']
    path = obtain(product['key'], product['size'], cache)
    with fits.open(path, memmap=False) as h:
        hdr, tab = h[1].header, h[1].data
        if h[0].header['TICID'] != tic or h[0].header['SECTOR'] != row['sector']:
            raise ValueError('archive metadata mismatch')
        cadence = round(86400 * hdr['TIMEDEL'])
        flag = tab['QUALITY'] == 0
        t = np.asarray(tab['TIME'][flag], float)
        f = np.asarray(tab['SAP_FLUX'][flag], float)
        ok = np.isfinite(t) & np.isfinite(f)
        t, f = t[ok], f[ok]
        result = search(t, f, cadence)
        match = None
        # Use literature only after the search. Require time drift across the visit
        # less than one full template width, and phase agreement at visit midpoint.
        # For the imprecise TIC 60040774 ephemeris use period alone in later sectors.
        period, epoch, _ = KNOWN[tic]
        for rank, peak in enumerate(result['periodic'], 1):
            p, e = peak['period_days'], peak['epoch_btjd']
            width = peak['duration_min']/1440
            middle = np.median(t)
            predicted = e + p * round((middle-e)/p)
            slip = abs((predicted-epoch+period/2) % period-period/2)
            uses_phase = tic != 60040774 or row['sector']==33
            if abs(p-period)*np.ptp(t)/p <= width and (slip <= width or not uses_phase):
                match = dict(rank=rank,phase_required=uses_phase,phase_offset_minutes=float(1440*slip),**peak)
                break
        return dict(tic=tic,sector=row['sector'],cadence_s=cadence,source_url=BASE+product['key'],
                    bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    columns_sha256=digest(tab),photometry='SAP_FLUX', crowding=float(hdr['CROWDSAP']), usable=result['usable_cadences'],
                    first_btjd=float(t.min()),last_btjd=float(t.max()),periodic_peaks=len(result['periodic']),
                    top=result['periodic'][0] if result['periodic'] else None,matching_peak=match)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inventory',type=Path)
    parser.add_argument('cache',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--previous-fetch',type=Path,
                        default=next((p/'results/tess-additional-transits/fetch.py' for p in HERE.parents
                                      if (p/'results/tess-additional-transits/fetch.py').exists()), None),
                        help='Path to repository results/tess-additional-transits/fetch.py')
    args=parser.parse_args()
    if args.previous_fetch is None:
        parser.error('Supply --previous-fetch with the repository helper path')
    digest = runpy.run_path(str(args.previous_fetch))['digest']
    args.cache.mkdir(parents=True, exist_ok=True)
    rows=json.loads(args.inventory.read_text())['rows']
    jobs=[(r, p) for r in rows for p in r['products'] if p['key'].endswith(('_lc.fits','_fast-lc.fits'))]
    def one(pair):
        row,product=pair
        try:
            return measure(row, product, args.cache, digest)
        except Exception as exc:
            return dict(tic=row['tic'],sector=row['sector'],source_url=BASE+product['key'],error=str(exc))
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        results=list(pool.map(one, jobs))
    args.output.write_text(json.dumps(results,indent=1,allow_nan=False)+'\n')
    print(f'{len(results)} products, {sum("error" in r for r in results)} errors')

if __name__ == '__main__':
    main()
