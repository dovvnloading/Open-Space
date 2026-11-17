"""Fetch the STScI catalogs in pairs.csv; cache only the columns used by the screen."""
import argparse
import csv
import hashlib
import time
from urllib.error import HTTPError, URLError
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
from pathlib import Path
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


def fetch(product, cache):
    url = product_url(product, 'cat.ecsv')
    for attempt in range(5):
        try:
            with urlopen(url, timeout=90) as response:
                data = response.read()
            break
        except (HTTPError, URLError, TimeoutError):
            if attempt == 4:
                raise
            time.sleep(2**attempt)
    # ECSV serializes SkyCoord as two degree columns. Read compact numerical columns.
    frame = pd.read_csv(BytesIO(data), sep=r'\s+', comment='#', usecols=COLUMNS)
    np.savez_compressed(cache / (product + '.npz'), **{k: frame[k].to_numpy() for k in COLUMNS})
    return [product, url, len(data), hashlib.sha256(data).hexdigest(), len(frame)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    with (HERE / 'pairs.csv').open() as stream:
        products = {row[k] for row in csv.DictReader(stream) for k in ('reference', 'search')}
    with ThreadPoolExecutor(5) as pool:
        receipts = list(pool.map(lambda p: fetch(p, args.cache), sorted(products)))
    with (HERE / 'catalog_inputs.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['product','url','bytes','sha256','sources'])
        writer.writerows(receipts)


if __name__ == '__main__':
    main()
