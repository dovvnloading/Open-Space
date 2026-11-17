"""Match dated unTimely catalogue entries to the propagated Phan templates."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'orbits'))
from propagate import epoch, family, sky

BASE = 'https://unwise.me/data/neo7/untimely-catalog/'
TILE = '0348m500'
INDEX = 'untimely_index-neo7.fits.gz'


def retrieve(url, path):
    if not path.exists():
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=60) as response:
                    raw = response.read()
                path.write_bytes(raw)
                break
            except OSError:
                if attempt == 2:
                    raise
    return hashlib.sha256(path.read_bytes()).hexdigest()


def acquire(cache):
    cache.mkdir(exist_ok=True, parents=True)
    folder = f'{BASE}{TILE[:3]}/{TILE}/'
    index_hash = retrieve(BASE + INDEX, cache / INDEX)
    page_hash = retrieve(folder, cache / 'listing.html')
    filenames = re.findall(r'href="([^"]+\.cat\.fits\.gz)"', (cache / 'listing.html').read_text())

    def one(name):
        digest = retrieve(folder + name, cache / name)
        return {'file': name, 'url': folder + name, 'sha256': digest}
    with ThreadPoolExecutor(8) as pool:
        rows = list(pool.map(one, filenames))
    (HERE / 'receipts.json').write_text(json.dumps({'index_url': BASE + INDEX,
        'index_sha256': index_hash, 'listing_sha256': page_hash, 'source_tile': TILE,
        'catalogues': rows}, indent=1) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('cache', type=Path)
    opt = p.parse_args()
    acquire(opt.cache)
