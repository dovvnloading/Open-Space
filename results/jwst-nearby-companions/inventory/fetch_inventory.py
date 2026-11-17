"""Fetch released JWST filesets and audit archive pagination; public data only."""
import argparse
import concurrent.futures
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import time
from pathlib import Path
from urllib.request import Request, urlopen

ENDPOINT = 'https://mast.stsci.edu/search/jwst/api/v0.1/search'
COLUMNS = ('fileSetName program title observtn visit targprop targname targ_ra targ_dec '
           'proposal_type exp_type opticalElements date_obs effexptm publicReleaseDate access '
           'is_psf bkgdtarg targcat targdesc obslabel').split()
SAMPLES = {
    'coron': [{'exp_type': 'NRC_CORON,MIR_4QPM,MIR_LYOT'}],
    'nearby_programs': [{'exp_type': 'MIR_IMAGE,NRC_IMAGE'}, {'program': '1184,1911,3964,4403,4857,6122,7564,7833,8142,8581,9079'}],
    'stellar_imaging': [{'exp_type': 'MIR_IMAGE,NRC_IMAGE'}, {'targcat': 'Star'}, {'effexptm': '>=300'}],
}

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def get_page(payload):
    req = Request(ENDPOINT, data=json.dumps(payload).encode(), headers={'Content-Type':'application/json'})
    for attempt in range(3):
        try:
            with urlopen(req, timeout=180) as response:
                body = response.read()
                return body, json.loads(body)
        except (OSError, ValueError) as error:
            if attempt == 2:
                raise
            print('Retry:', error, flush=True)
            time.sleep(2 ** (attempt + 1))

def run_query(name, options):
    out = options.out
    out.mkdir(parents=True, exist_ok=True)
    conditions = SAMPLES[name] + [{'publicReleaseDate': '<=' + options.as_of}]
    rows, pages = [], []
    seen = set()
    offset = 0
    while True:
        payload = {'conditions': conditions, 'limit': 1000, 'offset': offset, 'select_cols': COLUMNS,
                   'sort_by': ['fileSetName','ArchiveFileID'], 'sort_desc': [False, False]}
        body, reply = get_page(payload)
        batch = reply['results']
        pages.append({'retrieved_utc': utc(), 'request': payload, 'sha256_response': hashlib.sha256(body).hexdigest(),
                      'total_results': reply['totalResults'], 'returned': len(batch)})
        for row in batch:
            key = row['ArchiveFileID']
            if key in seen:
                raise ValueError(f'{name}: duplicate ArchiveFileID {key} across archive pages')
            seen.add(key)
        rows.extend(batch)
        offset += len(batch)
        print(name, offset, '/', reply['totalResults'], flush=True)
        if offset >= reply['totalResults']:
            break
        if not batch:
            raise ValueError(f'{name}: empty page before end of query')
    if any(p['total_results'] != len(rows) for p in pages):
        raise ValueError(f'{name}: result count changed during pagination')
    public = [row for row in rows if row['access'] == 'PUBLIC' and row['publicReleaseDate'].replace('T',' ') <= options.as_of]
    public.sort(key=lambda row: (row['publicReleaseDate'], row['date_obs'], row['fileSetName']), reverse=True)
    buf = io.StringIO(newline='')
    w = csv.DictWriter(buf, fieldnames=['ArchiveFileID'] + COLUMNS, extrasaction='ignore',lineterminator='\n')
    w.writeheader(); w.writerows(public)
    data = buf.getvalue().encode()
    (out / (name + '.csv.gz')).write_bytes(gzip.compress(data, mtime=0))
    (out / (name + '-receipt.json')).write_text(json.dumps({'endpoint':ENDPOINT,'meaning':'MAST JWST mission filesets; not detector-level FITS products',
        'returned':len(rows), 'public_at_cutoff':len(public),'cutoff_utc':options.as_of,'pages':pages,'sha256_csv':hashlib.sha256(data).hexdigest()},indent=2)+'\n')

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--as-of', default='2026-11-17 06:51:31', help='inclusive public release cutoff (UTC)')
    ap.add_argument('--out', type=Path, default=Path(__file__).parent / 'archive')
    ap.add_argument('samples', nargs='*', help='sample keys: ' + ', '.join(SAMPLES))
    options = ap.parse_args()
    samples = options.samples or list(SAMPLES)
    if any(name not in SAMPLES for name in samples):
        ap.error('unknown sample key')
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        list(executor.map(lambda name: run_query(name, options), samples))
