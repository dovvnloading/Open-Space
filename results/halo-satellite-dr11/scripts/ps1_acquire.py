"""Snapshot independent Pan-STARRS1 DR2 primary stack photometry in fixed cells."""
import concurrent.futures
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import time
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'comparison'
BASE = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/stack.csv'
COLS = 'objID,raStack,decStack,gPSFMag,rPSFMag,iPSFMag,iKronMag,gPSFMagErr,rPSFMagErr'
FIELDS = {'bootes': (209, 13.5), 'coma': (185.75, 22.9)}
NAMES = ['id', 'ra', 'dec', 'g', 'r', 'gerr', 'rerr']


def tile(item):
    field, i, j = item
    x0, y0 = FIELDS[field]
    lo, bot = x0 + i/2, y0 + j/2
    ra, dec = lo + .25, bot + .25
    # 0.355 deg encloses the corners even at the celestial equator.
    query = dict(ra=ra, dec=dec, radius=.355, primaryDetection=1,
                 columns=f'[{COLS}]', pagesize=100000)
    query.update({'gPSFMag.gte': 17, 'gPSFMag.lt': 22.5,
                  'rPSFMag.gte': 16, 'rPSFMag.lt': 22.5})
    url = BASE + '?' + urllib.parse.urlencode(query)
    receipt = dict(field=field, i=i, j=j, url=url)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=100) as response:
                raw = response.read()
            data = csv.DictReader(io.StringIO(raw.decode()))
            if data.fieldnames != COLS.split(','):
                raise ValueError(f'unexpected header {str(data.fieldnames)[:100]}')
            stars = {}
            count = 0
            for row in data:
                count += 1
                x,y,g,r,ipsf,ikron,eg,er = [float(row[c]) for c in COLS.split(',')[1:]]
                if not (lo <= x < lo + .5 and bot <= y < bot + .5 and 17 <= g < 22.5
                        and 16 <= r < 22.5 and 0 < eg < .15 and 0 < er < .15
                        and min(ipsf,ikron)>0 and abs(ipsf-ikron)<.05):
                    continue
                ident = row['objID']
                entry = [ident,x,y,g,r,eg,er]
                if ident not in stars or eg**2+er**2 < sum(e**2 for e in stars[ident][-2:]):
                    stars[ident] = entry
            if count >= query['pagesize']:
                raise ValueError('page limit reached')
            receipt.update(utc=dt.datetime.now(dt.timezone.utc).isoformat(), count=count,
                           stars=len(stars), raw_sha256=hashlib.sha256(raw).hexdigest())
            return receipt, list(stars.values())
        except Exception as exc:
            receipt['last_error'] = str(exc)
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(receipt)


def main():
    items = [(f,i,j) for f in FIELDS for i in range(6) for j in range(4)]
    OUT.mkdir(exist_ok=True)
    receipts, rows = [], {field: {} for field in FIELDS}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for receipt, data in pool.map(tile,items):
            receipts.append(receipt)
            for row in data:
                previous = rows[receipt['field']].get(row[0])
                if previous is None or sum(e**2 for e in row[-2:]) < sum(e**2 for e in previous[-2:]):
                    rows[receipt['field']][row[0]] = row
    for field, data in rows.items():
        text = io.StringIO(newline='')
        writer = csv.writer(text, lineterminator='\n')
        writer.writerow(NAMES)
        writer.writerows(data[k] for k in sorted(data))
        (OUT / f'ps1_{field}.csv.gz').write_bytes(gzip.compress(text.getvalue().encode(),mtime=0))
    (OUT/'ps1_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')


if __name__ == '__main__':
    main()
