"""Inventory SPOC time-series products for three externally chosen controls."""
import concurrent.futures as cf
import datetime as dt
import hashlib
import json
from pathlib import Path
import urllib.parse as p
import urllib.request as u
import xml.etree.ElementTree as ET

OUT = Path(__file__).resolve().parent
TICS = [267574918, 60040774, 219244444]
CUTOFF = 110  # Current sector on 2026-11-17, including partial deliveries.
BASE = 'https://stpubdata.s3.amazonaws.com/'
def one(job):
    tic, sector = job
    n = f'{tic:016d}'
    prefix = f'tess/public/tid/s{sector:04d}/' + '/'.join(n[i:i+4] for i in range(0,16,4)) + '/'
    url = BASE+'?'+p.urlencode({'list-type':2,'prefix':prefix})
    for retry in range(3):
        try:
            with u.urlopen(url, timeout=25) as response:
                raw = response.read()
            root = ET.fromstring(raw)
            items = [dict(key=row.findtext('{*}Key'), size=int(row.findtext('{*}Size'))) for row in root.findall('{*}Contents')]
            if root.findtext('{*}IsTruncated') != 'false':
                raise ValueError('truncated listing')
            return dict(tic=tic,sector=sector,url=url,status='ok',sha256=hashlib.sha256(raw).hexdigest(), products=[r for r in items if r['key'].endswith(('_lc.fits','_tp.fits','_fast-lc.fits','_fast-tp.fits'))])
        except Exception as exc:
            err=str(exc)
    return dict(tic=tic,sector=sector,url=url,status=err,products=[])
if __name__ == '__main__':
    jobs=[(tic,s) for s in range(CUTOFF,0,-1) for tic in TICS]
    with cf.ThreadPoolExecutor(max_workers=10) as e:
        results=list(e.map(one,jobs))
    (OUT/'control_inventory.json').write_text(json.dumps(dict(asof=dt.datetime.now(dt.timezone.utc).isoformat(),rows=results),indent=1)+'\n')
    for t in TICS:
        print(t,[(r['sector'],len(r['products'])) for r in results if r['tic']==t and r['products']],flush=True)
