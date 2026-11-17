"""Rebuild the reference-object photometry and compare the frozen magnitude cuts."""
import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import pathlib
import re
import urllib.request
from datetime import UTC, datetime

HERE = pathlib.Path(__file__).resolve().parent
SNR_FACTOR = 2.5 / math.log(10)
ZTF_IDS = {'AT2026nik': 'ZTF26aaxutnf', 'AT2026dbl': 'ZTF26aaglqak',
           'AT2026ahqq': 'ZTF26aadfszg'}
SOURCES = {name: f'https://api.alerce.online/ztf/v1/objects/{oid}/detections'
           for name, oid in ZTF_IDS.items()}
SOURCES['GRB261113A'] = 'https://gcn.nasa.gov/circulars/46263.json'
FIELDS = ['object', 'instrument', 'filter', 'mjd', 'mag', 'sigma', 'record']


def photometry(name, data):
    if name in ZTF_IDS:
        # Multiple candidate IDs at the same visit occur at CCD overlaps.
        # Retain the lowest photometric error per same time and band.
        unique = {}
        for row in data:
            if row['isdiffpos'] != 1 or not (0 < row['sigmapsf'] <= SNR_FACTOR / 5):
                continue
            key = (row['mjd'], row['fid'])
            if key not in unique or row['sigmapsf'] < unique[key]['sigmapsf']:
                unique[key] = row
        for (mjd, fid), row in sorted(unique.items()):
            yield [name, 'ZTF', str(fid), mjd, row['magpsf'], row['sigmapsf'], str(row['candid'])]
    else:
        # Use only this instrument and the circular's midpoints, not its exposure starts.
        trigger = datetime.fromisoformat('2026-11-13T14:31:02+00:00')
        mjd_trigger = (trigger - datetime(1858, 11, 17, tzinfo=UTC)).total_seconds() / 86400
        for line in data['body'].splitlines():
            match = re.fullmatch(r'(\S+)\s+(\d+\.\d+)\s+\S+\s+R\s+(\d+\.\d+)\s+\+/-\s+(\d+\.\d+)\s+\S+\s+AZT-33IK\s*', line)
            if match:
                start, hours, mag, sigma = match.groups()
                yield [name, 'AZT-33IK', 'R', mjd_trigger + float(hours) / 24, float(mag), float(sigma), start]


def select(rows):
    pairs = []
    for left, right in itertools.combinations(rows, 2):
        dt = float(right['mjd']) - float(left['mjd'])
        if left['filter'] != right['filter'] or not 0.02 <= dt <= 1.0:
            continue
        sig1, sig2 = float(left['sigma']), float(right['sigma'])
        if min(SNR_FACTOR / sig1, SNR_FACTOR / sig2) < 5:
            continue
        change = float(right['mag']) - float(left['mag'])
        significance = abs(change) / math.hypot(sig1, sig2)
        passes = abs(change) >= .30 and abs(change / dt) > 1 and significance >= 3
        pairs.append((passes, left['record'], right['record'], dt, change / dt, significance))
    return pairs


def rebuild():
    receipts = {}
    out = io.StringIO(newline='')
    writer = csv.writer(out, lineterminator='\n')
    writer.writerow(FIELDS)
    for name, url in SOURCES.items():
        with urllib.request.urlopen(url, timeout=40) as response:
            raw = response.read()
            receipts[name] = {'url': response.geturl(), 'utc': datetime.now(UTC).isoformat(),
                              'http_status': response.status, 'bytes': len(raw),
                              'sha256_http_body': hashlib.sha256(raw).hexdigest()}
        writer.writerows(photometry(name, json.loads(raw)))
    canonical = out.getvalue().encode('utf-8')
    (HERE / 'photometry.csv').write_bytes(canonical)
    receipts['photometry_csv_uncompressed_sha256'] = hashlib.sha256(canonical).hexdigest()
    (HERE / 'query_receipts.json').write_text(json.dumps(receipts, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch', action='store_true', help='Refresh inputs and query receipts from live services')
    if parser.parse_args().fetch:
        rebuild()
    with (HERE / 'photometry.csv').open() as file:
        rows = list(csv.DictReader(file))
    for name in SOURCES:
        obj = [row for row in rows if row['object'] == name]
        eligible = select(obj)
        passing = [p for p in eligible if p[0]]
        print(name, f'{len(obj)} epochs, {len(eligible)} pairs, {len(passing)} pass')
        if passing:
            best = max(passing, key=lambda p: p[-1])
            print('  strongest pair by significance', best[1:])
