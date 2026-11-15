"""Download the exact Chandra level-2 public event lists recorded in the manifest."""
import argparse
import csv
import hashlib
from pathlib import Path
from shutil import copyfileobj
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent


def fetch_products(cache):
    cache.mkdir(exist_ok=True, parents=True)
    for row in csv.DictReader((ROOT / 'chandra-products.csv').read_text().splitlines()):
        out = cache / (row['obsid'] + '.fits.gz')
        if out.exists():
            with out.open('rb') as handle:
                if hashlib.file_digest(handle, 'sha256').hexdigest() == row['sha256']:
                    continue
        tmp = out.with_suffix('.partial')
        with urlopen(row['url'], timeout=120) as src, tmp.open('wb') as dest:
            copyfileobj(src, dest)
        with tmp.open('rb') as handle:
            if hashlib.file_digest(handle, 'sha256').hexdigest() != row['sha256']:
                raise ValueError(f"archive digest differs for {row['obsid']}")
        tmp.replace(out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    fetch_products(parser.parse_args().cache)
