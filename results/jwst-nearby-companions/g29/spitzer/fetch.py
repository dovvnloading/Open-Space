"""Download four Spitzer SHA mosaics covering the G29-38 source and record their hashes."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
ENDPOINT = 'https://irsa.ipac.caltech.edu/SIA'
PARAMS = dict(POS='CIRCLE 352.193797 5.246759 0.004', BAND='0.0000034 0.000026',
              COLLECTION='spitzer_sha', FACILITY='Spitzer Space Telescope', RESPONSEFORMAT='CSV')
SELECTED = {('11124224', 'IRAC2'), ('11124224', 'IRAC4'), ('35008512', 'IRAC2'), ('10149376', 'MIPS24')}

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def download(row):
    url = row['access_url']
    output = HERE / 'cache' / url.rsplit('/', 1)[1]
    if not output.exists():
        raw = urlopen(url, timeout=150).read()
        if not raw.startswith(b'SIMPLE  ='):
            raise ValueError(url)
        output.write_bytes(raw)
    data = output.read_bytes()
    return {key: row[key] for key in ('obs_id', 'energy_bandpassname', 'access_url', 'target_name', 'proposal_id', 't_min', 't_max')} | {'bytes': len(data), 'sha256': sha(data)}

def fetch():
    (HERE / 'cache').mkdir(exist_ok=True)
    query_url = ENDPOINT + '?' + urlencode(PARAMS)
    response = urlopen(query_url, timeout=150).read()
    rows = list(csv.DictReader(io.StringIO(response.decode())))
    selected = [r for r in rows if r['access_url'].endswith('_maic.fits') and any(
        '/r' + aid + '/' in r['access_url'] and r['energy_bandpassname'] == band for aid, band in SELECTED)]
    if len(selected) != len(SELECTED):
        raise ValueError('Unexpected number of SHA science mosaics')
    with ThreadPoolExecutor(max_workers=4) as pool:
        inputs = list(pool.map(download, selected))
    out = dict(utc=datetime.now(timezone.utc).isoformat(), query_url=query_url, response_sha256=sha(response),
               response_rows=len(rows), files=inputs)
    (HERE / 'spitzer-query.json').write_text(json.dumps(out, indent=2) + '\n')
if __name__ == '__main__':
    fetch()
