"""Verify all frozen catalog bytes, and each cached FITS cube if available."""
import argparse
import csv
import hashlib
from pathlib import Path
root=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--require-cubes',action='store_true');args=p.parse_args()
with (root/'input_hashes.csv').open() as manifest:
    for item in csv.DictReader(manifest):
        folder=root if item['field']=='nep' else root/'wide'
        source=folder/item['file']
        if item['name']=='cube' and not source.exists() and not args.require_cubes:
            print('uncached',source);continue
        with source.open('rb') as f:
            digest=hashlib.file_digest(f,'sha256').hexdigest()
        if source.stat().st_size!=int(item['bytes']) or digest!=item['sha256']:
            raise ValueError(f"Input changed: {source}")
        print('verified',item['field'],item['name'],digest)
