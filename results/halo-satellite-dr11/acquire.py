"""Freeze Legacy Surveys DR11 southern Tractor photometry. Requires numpy, fitsio and Linux memfd."""
import argparse
import csv
import datetime
import gzip
import hashlib
import io
import json
import os
import shutil
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import fitsio
import numpy as np

ROOT = Path(__file__).resolve().parent
CAT = ROOT / 'catalog'
BASE = 'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr11'
INDEX = BASE + '/survey-bricks.fits.gz'
# BRIGHT, SATUR_G/R, ALLMASK_G/R, BAILOUT, MEDIUM, SGA large galaxy.
# Do not veto globular clusters (13) or nearby resolved dwarfs (17).
BAD_MASK = sum(1 << bit for bit in [1, 2, 3, 5, 6, 10, 11, 12])
COLUMNS = ('ls_id_dr11 release brickid objid brick_primary type maskbits fitbits ra dec '
           'flux_g flux_r flux_ivar_g flux_ivar_r mw_transmission_g mw_transmission_r '
           'nobs_g nobs_r fracflux_g fracflux_r fracmasked_g fracmasked_r fracin_g fracin_r '
           'psfdepth_g psfdepth_r ref_cat ref_id pmra pmdec pmra_ivar pmdec_ivar').split()
OUTPUT = ('id ra dec g r gerr rerr star clean type maskbits fitbits gdepth rdepth release brickid '
          'objid nobs_g nobs_r ref_cat ref_id pmra pmdec pmra_ivar pmdec_ivar').split()


def read_release_table(url, cols):
    """Stream into an anonymous in-memory file with a SHA-256 over unaltered HTTP bytes."""
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=100) as response, os.fdopen(os.memfd_create('dr11'), 'w+b') as f:
                if response.status != 200:
                    raise RuntimeError(f'Expected a complete HTTP 200 response; got {response.status}')
                digest = hashlib.sha256()
                while block := response.read(1024 * 1024):
                    digest.update(block)
                    f.write(block)
                size = f.tell()
                if response.headers.get('Content-Length') and int(response.headers['Content-Length']) != size:
                    raise ValueError(f'incomplete input {url}')
                f.flush()
                with fitsio.FITS(f'/proc/self/fd/{f.fileno()}') as hdus:
                    data = hdus[1].read(columns=cols)
                receipt = {'url': url, 'bytes': size, 'sha256': digest.hexdigest(),
                           'raw_rows': len(data), 'last_modified': response.headers.get('Last-Modified'),
                           'retrieved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
                return data, receipt
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def select(data, box):
    mag, err, depth = {}, {}, {}
    for band in ('g', 'r'):
        flux = data['flux_' + band].astype(float)
        ivar = data['flux_ivar_' + band].astype(float)
        transmit = data['mw_transmission_' + band].astype(float)
        with np.errstate(invalid='ignore', divide='ignore'):
            mag[band] = 22.5 - 2.5 * np.log10(flux / transmit)
            err[band] = 2.5 / np.log(10) / (flux * np.sqrt(ivar))
            depth[band] = 22.5 - 2.5 * np.log10(5 / np.sqrt(data['psfdepth_' + band]) / transmit)
    keep = ((data['ra'] >= box['ra'][0]) & (data['ra'] < box['ra'][1])
            & (data['dec'] >= box['dec'][0]) & (data['dec'] < box['dec'][1]) & data['brick_primary']
            & (mag['g'] >= 17) & (mag['g'] < 24) & np.isfinite(mag['r'])
            & np.isfinite(err['g']) & np.isfinite(err['r'])
            & (err['g'] > 0) & (err['r'] > 0))
    typ = np.char.strip(data['type'].astype(str))
    clean = (data['maskbits'] & BAD_MASK) == 0
    for band in ('g', 'r'):
        clean &= ((data['nobs_' + band] >= 1) & (err[band] <= 2.5 / np.log(10) / 5)
                  & (data['fracflux_' + band] < 5) & (data['fracmasked_' + band] < .4)
                  & (data['fracin_' + band] > .3))
    fields = [data['ls_id_dr11'], data['ra'], data['dec'], mag['g'], mag['r'], err['g'], err['r'],
              (typ == 'PSF').astype(int), clean.astype(int), typ, data['maskbits'], data['fitbits'],
              depth['g'], depth['r']]
    fields += [data[n] for n in OUTPUT[14:]]
    for row in zip(*(arr[keep] for arr in fields)):
        yield [f'{v:.7f}' if i in (1, 2) else f'{v:.5f}' if isinstance(v, np.floating)
               else str(v).strip() for i, v in enumerate(row)]


def process_brick(task):
    field, brick, box, temp = task
    result = temp / (brick + '.json')
    piece = temp / (brick + '.gz')
    if result.exists() and piece.exists():
        return json.loads(result.read_text())
    url = f'{BASE}/south/tractor/{brick[:3]}/tractor-{brick}.fits'
    data, receipt = read_release_table(url, COLUMNS)
    stream = io.StringIO()
    writer = csv.writer(stream, lineterminator='\n')
    count = 0
    for row in select(data, box):
        writer.writerow(row)
        count += 1
    payload = gzip.compress(stream.getvalue().encode(), compresslevel=3, mtime=0)
    piece.write_bytes(payload)
    receipt.update(field=field, brick=brick, selected_rows=count)
    result.write_text(json.dumps(receipt))
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', required=True, type=Path, help='Scratch space for resumeable compressed brick extracts')
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    fields = json.loads((CAT / 'footprint.json').read_text())['fields']
    index, receipt = read_release_table(INDEX, ['BRICKNAME', 'RA1', 'RA2', 'DEC1', 'DEC2'])
    (CAT / 'brick_index_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    tasks = []
    for field, box in fields.items():
        intersects = ((index['RA2'] > box['ra'][0]) & (index['RA1'] < box['ra'][1])
                      & (index['DEC2'] > box['dec'][0]) & (index['DEC1'] < box['dec'][1]))
        temp = args.work / field
        temp.mkdir(parents=True, exist_ok=True)
        for brick in sorted(index['BRICKNAME'][intersects]):
            tasks.append((field, str(brick).strip(), box, temp))
    print(f'{len(tasks)} intersecting bricks', flush=True)
    with ThreadPoolExecutor(args.workers) as pool:
        receipts = []
        for i, item in enumerate(pool.map(process_brick, tasks), 1):
            receipts.append(item)
            if i % 20 == 0:
                print(f'{i}/{len(tasks)} bricks', flush=True)
    with (CAT / 'brick_receipts.csv').open('w') as out:
        w = csv.DictWriter(out, fieldnames=list(receipts[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(receipts)
    summary = {}
    for field in fields:
        dst = CAT / (field + '.csv.gz')
        # Concatenated gzip members are RFC 1952 compliant; readers transparently concatenate.
        with dst.open('wb') as out:
            out.write(gzip.compress((','.join(OUTPUT) + '\n').encode(), mtime=0))
            for r in receipts:
                if r['field'] == field:
                    with (args.work / field / (r['brick'] + '.gz')).open('rb') as src:
                        shutil.copyfileobj(src, out)
        expected = sum(r['selected_rows'] for r in receipts if r['field'] == field)
        ids = set()
        count = stars = 0
        with gzip.open(dst, 'rt') as f:
            for row in csv.DictReader(f):
                ids.add(int(row['id']))
                count += 1
                stars += int(row['star']) * int(row['clean'])
        if expected != count or len(ids) != count:
            raise ValueError(f'{field}: expected {expected} rows, got {count} with {len(ids)} unique IDs')
        summary[field] = {'rows': count, 'unique_ids': len(ids), 'clean_stars': stars,
                          'bricks': sum(r['field'] == field for r in receipts),
                          'sha256': hashlib.sha256(dst.read_bytes()).hexdigest()}
    (CAT / 'snapshot.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
