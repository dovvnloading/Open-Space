"""Fetch the frozen MAST products and verify time-series columns used in the screen."""
import argparse
import csv
import hashlib
import io
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
COLUMNS = ('TIME', 'PDCSAP_FLUX', 'PDCSAP_FLUX_ERR', 'SAP_FLUX', 'SAP_BKG',
           'QUALITY', 'MOM_CENTR1', 'MOM_CENTR2')


def digest(columns):
    result = hashlib.sha256()
    for name in COLUMNS:
        data = np.asarray(columns[name], dtype='<i8' if name == 'QUALITY' else '<f8').copy()
        if name != 'QUALITY':
            data[np.isnan(data)] = np.nan
        result.update(data.tobytes())
    return result.hexdigest()


def fetch(product, cache):
    path = cache / f"{product['tic']}-{product['sector']}.npz"
    expected = product['columns_sha256']
    if path.exists():
        with np.load(path) as saved:
            if digest(saved) == expected:
                return
    with urlopen(product['url'], timeout=120) as response:
        with fits.open(io.BytesIO(response.read()), memmap=False) as hdus:
            columns = {name: np.asarray(hdus[1].data[name]) for name in COLUMNS}
            if digest(columns) != expected:
                raise ValueError(f"Archive columns changed for {product['url']}")
            temporary = path.with_suffix('.pending.npz')
            np.savez_compressed(temporary, **columns)
            temporary.replace(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    with (HERE / 'products.csv').open() as source:
        products = list(csv.DictReader(source))
    with ThreadPoolExecutor(max_workers=6) as pool:
        list(pool.map(lambda row: fetch(row, args.cache), products))
    print(f"Verified {len(products)} MAST products")
