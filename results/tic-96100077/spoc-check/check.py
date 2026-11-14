"""Repeat the target-pixel clock measurement on two published SPOC FFI light curves."""
import argparse
import csv
import hashlib
import io
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io import fits

HERE = Path(__file__).resolve().parent
COLUMNS = ('TIME', 'PDCSAP_FLUX', 'QUALITY')


def column_digest(data):
    h = hashlib.sha256()
    for key in COLUMNS:
        arr = np.array(data[key], dtype='<i8' if key == 'QUALITY' else '<f8', copy=True)
        if key != 'QUALITY':
            arr[np.isnan(arr)] = np.nan
        h.update(arr.tobytes())
    return h.hexdigest()


def sector_light_curve(row, cache):
    sector = int(row['sector'])
    filename = cache / f'TESS-SPOC-{sector}.fits'
    if not filename.exists():
        with urlopen(row['url'], timeout=90) as response:
            payload = response.read()
        if hashlib.sha256(payload).hexdigest() != row['fits_sha256']:
            raise ValueError(f'Archive FITS changed for sector {sector}')
        pending = filename.with_suffix('.pending')
        pending.write_bytes(payload)
        pending.replace(filename)
    payload = filename.read_bytes()
    if hashlib.sha256(payload).hexdigest() != row['fits_sha256']:
        raise ValueError(f'FITS digest mismatch for sector {sector}')
    with fits.open(io.BytesIO(payload), memmap=False) as hdus:
        data = hdus[1].data
        if column_digest(data) != row['columns_sha256']:
            raise ValueError(f'Column digest mismatch for sector {sector}')
        good = (data['QUALITY'] == 0) & np.isfinite(data['TIME']) & np.isfinite(data['PDCSAP_FLUX'])
        flux = data['PDCSAP_FLUX'][good].astype(float)
        return pd.DataFrame(dict(sector=sector, time=data['TIME'][good], flux=flux/np.median(flux)))


if __name__ == '__main__':
    import sys
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--analysis-dir', type=Path, default=HERE.parent,
                        help='Directory containing the historical clock measurement code, measure.py')
    args = parser.parse_args()
    sys.path.insert(0, str(args.analysis_dir.resolve()))
    from measure import windows, fit
    args.cache.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    with (HERE / 'products.csv').open() as file:
        curves = pd.concat([sector_light_curve(row, args.cache) for row in csv.DictReader(file)])
    timings = windows(curves)
    timings.to_csv(args.output / 'timings.csv', index=False, float_format='%.10f')
    rows = []
    for label in ('both', 'primary', 'secondary'):
        for local in (False, True):
            _, baseline, _ = fit(timings, method=label, use_orbit=False, local=local)
            coeff, fitted, _ = fit(timings, method=label, use_orbit=True, local=local)
            rows.append(dict(clock=label, sector_offsets=local,
                             n_windows=int(timings[label].notna().sum()),
                             baseline_rms_s=baseline, orbit_rms_s=fitted, orbit_scale=coeff[-1]))
    pd.DataFrame(rows).to_csv(args.output / 'summary.csv', index=False, float_format='%.8f')
