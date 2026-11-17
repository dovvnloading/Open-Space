"""Download ALeRCE's ZTF discovery cohort, independently counting rows per UTC day.

Run with Python 3.10+ from any directory. The query includes negative detections
so that signed differences, ignored in slope measurements, remain auditable.
"""
import concurrent.futures
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import time
import urllib.error
import urllib.parse
import urllib.request

ROOT = Path(os.environ['FAST_ALERTS_NEW_CAMPAIGN']).resolve() if os.environ.get('FAST_ALERTS_NEW_CAMPAIGN') else Path(__file__).resolve().parents[1]
RAW = ROOT / 'data' / 'ztf'
META = ROOT / 'receipts' / 'ztf'
TAP = 'https://tap.alerce.online/tap/sync'
WINDOW = (61345, 61361)
COLUMNS = '''o.oid,o.firstmjd,o.lastmjd,o.ndet,o.mjdstarthist,o.stellar,
o.meanra,o.meandec,d.candid,d.pid,d.mjd,d.fid,d.magpsf,d.sigmapsf,
d.isdiffpos,d.drb,d.rb,d.ra,d.dec,d.has_stamp,d.dubious,d.rfid,d.diffmaglim'''
# CSV avoids loss of 64-bit integer candidate/processing IDs.
def get(tag, q, bodydir=RAW):
    if not os.environ.get('FAST_ALERTS_NEW_CAMPAIGN') or Path(os.environ['FAST_ALERTS_NEW_CAMPAIGN']).resolve() != ROOT.resolve():
        raise RuntimeError('Live retrieval requires a separate exclusive FAST_ALERTS_NEW_CAMPAIGN directory')
    if ROOT.is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('Live responses must be outside the frozen result')
    if (bodydir / (tag + '.csv.gz')).exists() or (META / (tag + '.json')).exists():
        raise FileExistsError(tag)
    params = dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', MAXREC='200000', QUERY=q)
    url = TAP + '?' + urllib.parse.urlencode(params)
    last = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=150) as response:
                raw = response.read()
                if not response.headers.get('Content-Type', '').startswith('text/csv'):
                    raise ValueError(str(raw[:1000]))
                info = dict(method='GET', url=url, query=q, http_status=response.status,
                            retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                            mime=response.headers.get('Content-Type'), bytes=len(raw),
                            sha256=hashlib.sha256(raw).hexdigest())
            dest = bodydir / (tag + '.csv.gz')
            dest.write_bytes(gzip.compress(raw, mtime=0))
            info['file'] = str(dest.relative_to(ROOT))
            (META / (tag + '.json')).write_text(json.dumps(info, indent=2) + '\n')
            return list(csv.DictReader(io.StringIO(raw.decode())))
        except (OSError, ValueError) as exc:
            last = exc
            time.sleep(1 + attempt)
    raise RuntimeError(f'{tag} failed: {last}')
def day(n):
    clause = (f' FROM ztf.object AS o INNER JOIN ztf.detection AS d ON o.oid = d.oid'
              f' WHERE o.firstmjd >= {n} AND o.firstmjd < {n+1} AND o.ndet >= 2'
              f' AND d.mjd >= {WINDOW[0]} AND d.mjd < {WINDOW[1]}')
    count = int(get(f'{n}_count', 'SELECT COUNT(*) AS nrows'+clause, META)[0]['nrows'])
    rows = get(f'{n}_detections','SELECT '+ COLUMNS.replace('\n',' ') + clause)
    if len(rows) != count:
        raise ValueError(f'{n}: {count} expected; received {len(rows)}')
    print(f'{n} {len(rows)}',flush=True)
def main():
    if not os.environ.get('FAST_ALERTS_NEW_CAMPAIGN'):
        raise RuntimeError('Set FAST_ALERTS_NEW_CAMPAIGN to a new empty directory')
    if ROOT.is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('New acquisition must be outside the frozen result')
    ROOT.mkdir(exist_ok=False, parents=True)
    RAW.mkdir(exist_ok=False, parents=True)
    META.mkdir(exist_ok=False, parents=True)
    get('columns', "SELECT table_name,column_name,description FROM TAP_SCHEMA.columns WHERE table_name IN ('ztf.object','ztf.detection','ztf.magstat','ztf.ps1_ztf','ztf.gaia_ztf','ztf.ss_ztf')", META)
    get('objects_total', 'SELECT COUNT(*) AS nobj FROM ztf.object WHERE firstmjd >= 61345 AND firstmjd < 61361 AND ndet >= 2', META)
    get('latest_object', 'SELECT TOP 1 firstmjd FROM ztf.object ORDER BY firstmjd DESC', META)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(day, range(*WINDOW)))
if __name__ == '__main__':
    main()
