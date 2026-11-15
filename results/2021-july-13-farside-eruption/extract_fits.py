"""Rebuild the count-rate and geometry extracts from the mission FITS products.

Requires numpy and astropy. Downloads to a temporary directory by default.
For an offline rerun: python extract_fits.py --cache /path/to/downloaded/fits
"""
import argparse
import csv
from datetime import datetime, timedelta
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlretrieve
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--cache', type=Path)
args = parser.parse_args()
manifest = json.loads((ROOT/'inputs/sources.json').read_text())
rates = []
geometry = {}
with TemporaryDirectory() as temporary:
    cache = args.cache or Path(temporary)
    cache.mkdir(parents=True, exist_ok=True)
    for name, meta in manifest.items():
        if not name.endswith('.fits'):
            continue
        source = cache / name
        if not source.exists():
            urlretrieve(meta['url'], source)
        if hashlib.sha256(source.read_bytes()).hexdigest() != meta['sha256']:
            raise ValueError(f'Checksum mismatch for {name}')
        with fits.open(source) as fitsfile:
            head = fitsfile[0].header
            data = fitsfile['DATA'].data
            bands = fitsfile['ENERGIES'].data
            date = datetime.fromisoformat(head['DATE-OBS'])
            if any(fitsfile['DATA'].columns[k].unit != 'cs' for k in ['time','timedel']):
                raise ValueError('Unexpected telemetry time unit')
            centers = [date + timedelta(seconds=t/100.) for t in data['time']]
            durations = data['timedel'] / 100.
            counts = data['counts'].reshape(len(data), -1, len(bands)).sum(axis=1)
            instrument = 'SCI' if 'sci-xray' in name else 'QL'
            geometry[instrument] = {k: head[k] for k in ['DATE-AVG','RSUN_ARC','HGLT_OBS','HGLN_OBS','CRLN_OBS','EAR_TDEL']}
            for low, high in [(4,10),(15,25)]:
                channels = (bands['e_low'] >= low) & (bands['e_high'] <= high)
                rate = counts[:,channels].sum(axis=1) / durations
                for tm, val in zip(centers, rate):
                    if datetime(2021,7,13,17,40) <= tm < datetime(2021,7,13,19,5):
                        rates.append([instrument,f'{low}-{high}',tm.isoformat(timespec='milliseconds'),f'{val:.6f}'])
with (ROOT/'inputs/stix_rates.csv').open('w') as out:
    writer = csv.writer(out)
    writer.writerow(['product','band_keV','spacecraft_utc','count_rate_s'])
    writer.writerows(sorted(rates))
(ROOT/'inputs/stix_headers.json').write_text(json.dumps(geometry,indent=2,sort_keys=True)+'\n')
