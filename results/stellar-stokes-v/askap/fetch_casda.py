"""Fetch the RACS-low2 release Stokes V detections, column definitions and observing epochs."""
import concurrent.futures
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
URL = 'https://casda.csiro.au/casda_vo_tools/tap/sync'
V_TABLE = 'AS110.racs_low2_v_v01'


def download(name, query):
    url = URL + '?' + urlencode({'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv',
                                 'MAXREC': 10000, 'QUERY': query})
    with urlopen(url, timeout=180) as response:
        data = response.read()
    if data.startswith(b'<'):
        raise ValueError(data[:800])
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    if not rows or len(rows) >= 10000:
        raise ValueError(f'Empty or truncated {name}: {len(rows)}')
    (ROOT / (name + '.csv')).write_bytes(data)
    return dict(name=name, url=url, query=query, rows=len(rows), bytes=len(data),
                sha256=hashlib.sha256(data).hexdigest(), utc=datetime.now(timezone.utc).isoformat())


if __name__ == '__main__':
    queries = {'racs_low2_v': f'SELECT * FROM {V_TABLE} ORDER BY id',
               'racs_low2_schema': 'SELECT column_name,description,unit FROM TAP_SCHEMA.columns '
                                   f"WHERE table_name='{V_TABLE}' ORDER BY column_name"}
    with concurrent.futures.ThreadPoolExecutor(2) as pool:
        receipts = list(pool.map(lambda item: download(*item), queries.items()))
    (ROOT / 'casda_receipts.json').write_text(json.dumps(receipts, indent=2) + '\n')
