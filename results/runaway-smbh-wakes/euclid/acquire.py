"""Freeze the public Q1 MER products and retrieve preselected native tile-centre cutouts.

Run from any directory. Output tables go under --out and FITS under --cache.
Only the compact receipts and tables need version control.
"""
import argparse
import csv
import hashlib
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

import numpy as np

S3 = 'https://nasa-irsa-euclid-q1.s3.us-east-1.amazonaws.com/'
IBE = 'https://irsa.ipac.caltech.edu/ibe/data/euclid/'
INDEX = S3 + '?' + urlencode({'list-type': 2, 'prefix': 'q1/MER/', 'delimiter': '/'})
CHECKSUMS = S3 + 'q1/MER/euclid_q1-MER.md5'
# Tiles fixed before looking at these VIS pixels, copied from the earlier Q1 search.
PRESELECTED = {'102020065', '102022484', '102044820', '102159484'}
TERMS = {'VIS': 'BGSUB-MOSAIC-VIS_', 'rms': 'MOSAIC-VIS-RMS_',
         'Y': 'BGSUB-MOSAIC-NIR-Y_', 'J': 'BGSUB-MOSAIC-NIR-J_', 'H': 'BGSUB-MOSAIC-NIR-H_'}
FITS_TYPES = {'-32': '>f4', '-64': '>f8', '16': '>i2', '32': '>i4', '8': 'u1'}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url, path, header=False):
    if path.is_file():
        return path.read_bytes()
    request = Request(url, headers={'Range': 'bytes=0-5759'} if header else {})
    for attempt in range(3):
        try:
            with urlopen(request, timeout=75) as response:
                if header and (response.status != 206 or
                               not response.headers.get('Content-Range', '').startswith('bytes 0-5759/')):
                    raise ValueError('unexpected archived FITS byte range')
                data = response.read()
                if header and len(data) != 5760:
                    raise ValueError('incomplete FITS header')
                if (header or path.suffix == '.fits') and data[:8] != b'SIMPLE  ':
                    raise ValueError('unexpected data type from archive')
            tmp = path.with_suffix(path.suffix + '.tmp')
            tmp.write_bytes(data)
            tmp.replace(path)
            return data
        except (OSError, ValueError):
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def parse_primary(raw, image=True):
    """Read the uncompressed survey primary HDU including FITS physical scaling."""
    h = {}
    for pos in range(0, len(raw), 80):
        card = raw[pos:pos+80].decode('ascii')
        key = card[:8].strip()
        if key == 'END':
            break
        if card[8:10] == '= ':
            value = card[10:].strip()
            h[key] = value.split("'")[1].strip() if value.startswith("'") else value.split('/')[0].strip()
    else:
        raise ValueError('FITS END is missing')
    if h.get('SIMPLE') != 'T' or h.get('NAXIS') != '2':
        raise ValueError('unexpected FITS primary header')
    if not image:
        return h
    shape = int(h['NAXIS2']), int(h['NAXIS1'])
    off = (pos + 80 + 2879) // 2880 * 2880
    a = np.frombuffer(raw, FITS_TYPES[h['BITPIX']], count=shape[0]*shape[1], offset=off).reshape(shape)
    return h, a.astype('f4') * float(h.get('BSCALE', 1)) + float(h.get('BZERO', 0))


def write_csv(path, rows):
    with path.open('w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def inventory(root, cache):
    listing = fetch(INDEX, cache / 'tiles.xml')
    e = ET.fromstring(listing)
    if e.findtext('{*}IsTruncated') != 'false':
        raise ValueError('S3 directory pagination required')
    tiles = sorted(x.text.split('/')[2] for x in e.findall('{*}CommonPrefixes/{*}Prefix'))
    checksums = fetch(CHECKSUMS, cache / 'mer.md5')
    if hashlib.md5(checksums).hexdigest() != e.findtext('{*}Contents/{*}ETag').strip('"'):
        raise ValueError('checksum index differs from directory ETag')
    products = {tile: {} for tile in tiles}
    for line in checksums.decode().splitlines():
        md5, rel = line.split()
        for band, term in TERMS.items():
            if term in rel:
                tile = rel.split('/')[0]
                if band in products[tile] or not re.fullmatch('[0-9a-f]{32}', md5):
                    raise ValueError(rel)
                products[tile][band] = (rel, md5)
    rows = []
    for tile, bands in products.items():
        if set(bands) != set(TERMS):
            raise ValueError((tile, list(bands)))
        for band, (rel, md5) in sorted(bands.items()):
            rows.append(dict(tile=tile, band_or_product=band, relative_to_mer=rel, publisher_md5=md5))
    write_csv(root / 'q1_products.csv', rows)
    (root / 'archive_queries.json').write_text(json.dumps([
        dict(url=INDEX, bytes=len(listing), sha256=sha(listing), tile_ids=len(tiles),
             is_truncated=False),
        dict(url=CHECKSUMS, bytes=len(checksums), sha256=sha(checksums),
             md5=hashlib.md5(checksums).hexdigest(), archive_lines=len(checksums.splitlines()),
             selected_product_rows=len(rows))], indent=2) + '\n')
    return products


def get_center(tile, bands, cache):
    url = S3 + 'q1/MER/' + bands['VIS'][0] + '?range_start=0&range_end=5759'
    raw = fetch(url, cache / (tile+'.head'), header=True)
    h = parse_primary(raw, False)
    if h.get('FILTER') != 'VIS' or min(int(h['NAXIS1']), int(h['NAXIS2'])) < 1025:
        raise ValueError('unexpected VIS tile header')
    ra, dec = float(h['CRVAL1']), float(h['CRVAL2'])
    # The four disconnected target groups are widely separated in declination.
    field = 'EDF-N' if dec > 55 else 'EDF-S' if dec < -40 else 'EDF-F' if dec < -20 else 'LDN1641'
    return dict(tile=tile, ra_deg=ra, dec_deg=dec, field=field, vis_magzero=float(h['MAGZERO']),
                header_bytes=len(raw), header_sha256=sha(raw))


def download(row, band, bands, cache):
    tile = row['tile']
    size = 102.4 if tile in PRESELECTED else 51.2
    url = (IBE + 'q1/MER/' + bands[band][0] + f"?center={row['ra_deg']:.8f},{row['dec_deg']:.8f}deg"
           + f'&size={size}arcsec&gzip=false')
    dest = cache / f'{tile}_{band}.fits'
    data = fetch(url, dest)
    h, a = parse_primary(data)
    if h['FILTER'] != ('VIS' if band == 'rms' else band) and h['FILTER'] != 'NIR_'+band:
        # Some NISP MER tiles record filters with a NISP_ prefix.
        if band not in ('Y', 'J', 'H') or not h['FILTER'].endswith(band):
            raise ValueError('unexpected band '+h['FILTER'])
    if a.shape[0] < size * 9.5 or a.shape[1] < size * 9.5:
        raise ValueError(f'truncated cutout {a.shape}')
    good = np.isfinite(a) & (a != 0)
    if band == 'rms':
        good &= (a > 0) & (a < 1e8)
    return dict(tile=tile, field=row['field'], product=band, size_arcsec=size, nx=a.shape[1], ny=a.shape[0],
                pixels_passing_single_product_mask=int(good.sum()), bytes=len(data), sha256=sha(data),
                url=url)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=6)
    parser.add_argument('--inventory-only', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)
    args.cache.mkdir(exist_ok=True, parents=True)
    products = inventory(args.out, args.cache)
    with ThreadPoolExecutor(args.jobs) as pool:
        centers = list(pool.map(lambda t: get_center(t, products[t], args.cache), products))
        write_csv(args.out / 'q1_tiles.csv', centers)
        if args.inventory_only:
            return
        tasks = [(r, b) for r in centers if r['field'] != 'LDN1641' for b in
                 (TERMS if r['tile'] in PRESELECTED else ('VIS', 'rms'))]
        def run(task):
            row, band = task
            try:
                receipt = download(row, band, products[row['tile']], args.cache)
                return receipt, None
            except (OSError, ValueError) as exc:
                return None, dict(tile=row['tile'], product=band, error=repr(exc))
        pairs = list(pool.map(run, tasks))
    write_csv(args.out / 'q1_cutouts.csv', [r for r, err in pairs if r])
    (args.out / 'q1_errors.json').write_text(json.dumps([err for _, err in pairs if err], indent=2)+'\n')


if __name__ == '__main__':
    main()
