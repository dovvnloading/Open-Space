"""Fetch and verify frozen STScI catalog bytes before preparing replay inputs."""
import argparse
import csv
import hashlib
import re
import time
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
COLUMNS = ['label', 'sky_centroid.ra', 'sky_centroid.dec', 'aper_total_flux',
           'aper_total_flux_err', 'is_extended', 'nn_dist', 'aper50_flux', 'aper70_flux']


def product_url(product, suffix):
    program, visit = product.split('_')[0].split('-')
    target_type = product.split('_')[1][0]
    return f'https://stpubdata.s3.amazonaws.com/jwst/public/{program}/L3/{target_type}/{visit}/{product}_{suffix}'


def frozen_inputs(directory=HERE):
    with (directory / 'pairs.csv').open(newline='', encoding='utf-8') as stream:
        products = {row[k] for row in csv.DictReader(stream) for k in ('reference', 'search')}
    with (directory / 'catalog_inputs.csv').open(newline='', encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    manifest = {}
    for row in rows:
        name = row['product']
        if (not name or Path(name).name != name or name in manifest or
                row['url'] != product_url(name, 'cat.ecsv') or
                not re.fullmatch(r'[0-9a-f]{64}', row['sha256']) or int(row['bytes']) < 1):
            raise ValueError(f'invalid frozen catalog receipt: {name!r}')
        manifest[name] = row
    if set(manifest) != products:
        raise ValueError('catalog receipt products do not exactly match the frozen pairs')
    return manifest


def verified_bytes(item, cache):
    name, row = item
    raw_file = cache / (name + '.cat.ecsv')
    if raw_file.exists():
        data = raw_file.read_bytes()
    else:
        for attempt in range(5):
            try:
                with urlopen(row['url'], timeout=90) as response:
                    data = response.read()
                break
            except (HTTPError, URLError, TimeoutError):
                if attempt == 4:
                    raise
                time.sleep(2**attempt)
    if len(data) != int(row['bytes']) or hashlib.sha256(data).hexdigest() != row['sha256']:
        raise ValueError(f'frozen catalog byte mismatch: {name}; receipt and outputs unchanged')
    return name, data, not raw_file.exists()


def prepare(cache, directory=HERE):
    """Validate all products, then parse; existing NPZ files are never trusted as input."""
    manifest = frozen_inputs(directory)
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=5) as pool:
        raw = list(pool.map(lambda item: verified_bytes(item, cache), sorted(manifest.items())))
    # No derived output is touched until every product is verified and parsed.
    parsed = []
    for name, data, missing in raw:
        frame = pd.read_csv(BytesIO(data), sep=r'\s+', comment='#', usecols=COLUMNS)
        if len(frame) != int(manifest[name]['sources']):
            raise ValueError(f'frozen catalog source count mismatch: {name}')
        parsed.append((name, data, missing, {k: frame[k].to_numpy() for k in COLUMNS}))
    for name, data, missing, columns in parsed:
        if missing:
            (cache / (name + '.cat.ecsv')).write_bytes(data)
        np.savez_compressed(cache / (name + '.npz'), **columns)
    print(f'verified {len(parsed)} frozen catalogs; receipt unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    args = parser.parse_args()
    prepare(args.cache)


if __name__ == '__main__':
    main()
