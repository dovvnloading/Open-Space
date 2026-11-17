"""Verify frozen control inputs before deriving Gaia and published RV tables."""
import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
import tempfile
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

HERE = Path(__file__).resolve().parent
ESA_SHA = 'c5c7165bd13ca7bdad0454719b10c5637b4b2d32'
PAPERS = {'BH1': ('2209.06833v3', 41), 'BH2': ('2302.07880v2', 48), 'NS1': ('2402.06722v2', 59)}
IDS = {'4373465352415301632', '5870569352746779008', '6328149636482597888'}
QUERY = ('SELECT * FROM gaiadr3.nss_two_body_orbit WHERE source_id IN '
         '(4373465352415301632,5870569352746779008,6328149636482597888)')
GAIA_URL = 'https://gea.esac.esa.int/tap-server/tap/sync?' + urlencode(
    dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=QUERY))
ROW = re.compile(r'^\s*(24\d{5}\.\d+)\s+([−\-\d.]+)\s*±\s*([\d.]+)\s+([\w-]+)')
EXPECTED_URLS = {'dr3.csv': GAIA_URL, **{name.lower()+'.pdf':f'https://arxiv.org/pdf/{paper}'
                 for name,(paper,_) in PAPERS.items()},
                 **{name:f'https://raw.githubusercontent.com/esa/gaia-bhthree/{ESA_SHA}/{name}'
                    for name in ['epochast.dat','epochrv.dat']}}


def frozen_receipts(path=HERE/'sources.json'):
    data = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if set(data) != set(EXPECTED_URLS):
        raise ValueError('frozen source membership mismatch')
    for name, url in EXPECTED_URLS.items():
        row = data[name]
        if (row.get('url') != url or not re.fullmatch(r'[0-9a-f]{64}', row.get('sha256',''))
                or not isinstance(row.get('bytes'),int) or row['bytes'] <= 0):
            raise ValueError(f'invalid frozen receipt: {name}')
    return data


def verified_inputs(cache, receipts):
    raw = {}
    for name,row in receipts.items():
        path = cache/name
        if path.exists():
            data = path.read_bytes()
        else:
            with urlopen(row['url'],timeout=90) as response:
                data = response.read()
        if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ValueError(f'frozen source byte mismatch: {name}; receipts and tables unchanged')
        raw[name] = data
    return raw


def parse_inputs(raw, scratch):
    reader = csv.DictReader(io.StringIO(raw['dr3.csv'].decode('utf-8-sig'), newline=''))
    records = sorted(reader, key=lambda r:r['source_id'])
    if len(records) != len(IDS) or {r['source_id'] for r in records} != IDS:
        raise ValueError('incomplete or duplicate Gaia control set')
    for name, min_columns in [('epochast.dat', 8), ('epochrv.dat', 4)]:
        array = np.loadtxt(io.BytesIO(raw[name]), ndmin=2)
        if array.shape[0] == 0 or array.shape[1] < min_columns or not np.isfinite(array[:, :min_columns]).all():
            raise ValueError(f'invalid Gaia epoch table: {name}')
    all_rows = []
    for name,(_,expected) in PAPERS.items():
        pdf = scratch/(name.lower()+'.pdf')
        pdf.write_bytes(raw[name.lower()+'.pdf'])
        lines = subprocess.check_output(['pdftotext','-layout',str(pdf),'-'],text=True).splitlines()
        rows = [[name,t,v.replace('−','-'),err,instrument] for line in lines
                if (m := ROW.match(line)) for t,v,err,instrument in [m.groups()]]
        if len(rows) != expected or len({r[1] for r in rows}) != expected:
            raise ValueError(f'unexpected RV row count for {name}: {len(rows)}')
        all_rows.extend(rows)
    return reader.fieldnames, records, all_rows


def fetch(cache, out=None, receipt=HERE/'sources.json'):
    cache = Path(cache)
    out = Path(out) if out is not None else cache/'prepared'
    if out.resolve() == HERE.resolve() or HERE.resolve() in out.resolve().parents:
        raise ValueError('output must be separate from the frozen source directory')
    receipts = frozen_receipts(receipt)  # Validate membership and URLs before reading data.
    cache.mkdir(parents=True,exist_ok=True)
    raw = verified_inputs(cache,receipts)  # Verify every byte before parsing any source.
    with tempfile.TemporaryDirectory(dir=cache) as tmp:
        fields,records,rows = parse_inputs(raw,Path(tmp))
    # Only after all inputs parse successfully may any cache or output be written.
    out.mkdir(parents=True,exist_ok=True)
    for name,data in raw.items():
        path = cache/name
        if not path.exists():path.write_bytes(data)
    with (out/'nss_controls.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(records)
    with (out/'published_rv.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream);writer.writerow(['control','hjd_utc','rv_kms','err_kms','instrument']);writer.writerows(rows)
    print(f'verified {len(raw)} frozen control inputs; receipt unchanged; tables in {out}')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache',type=Path,help='Directory for verified source bytes')
    parser.add_argument('--out',type=Path,help='Separate table output directory (default: CACHE/prepared)')
    a=parser.parse_args();fetch(a.cache,a.out)
