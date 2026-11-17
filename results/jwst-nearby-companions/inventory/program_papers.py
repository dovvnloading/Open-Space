"""Get MAST's curated bibliography for each public program in visits.csv."""
import argparse
import concurrent.futures
import csv
import datetime as dt
import hashlib
import io
import json
from pathlib import Path
import time
from urllib.request import Request, urlopen

ENDPOINT='https://mast.stsci.edu/papertrack/api/v0.1/public-search'
FIELDS=['program','bibcode','first_author','title','publication_date','doi','papertypes']

def papers(program):
    params={'pagesize':10000,'mission':'JWST','jwst_programs':str(program),'format':'csv'}
    req=Request(ENDPOINT,data=json.dumps(params).encode(),headers={'Content-Type':'application/json'})
    for attempt in range(3):
        try:
            with urlopen(req,timeout=50) as response:
                data=response.read()
            break
        except OSError:
            if attempt==2:
                raise
            time.sleep(2**attempt)
    reader=csv.DictReader(io.StringIO(data.decode('utf-8-sig')))
    # Empty results are emitted as an empty body or CSV header by Papertrack.
    if data.strip() and (not reader.fieldnames or 'bibcode' not in reader.fieldnames):
        raise ValueError(f'Unexpected response for {program}: {data[:150]!r}')
    rows=list(reader)
    if len(rows)>=params['pagesize']:
        raise ValueError('Bibliography page limit reached')
    receipt={'query':params,'sha256_response':hashlib.sha256(data).hexdigest(),'rows':len(rows),
             'retrieved_utc':dt.datetime.now(dt.timezone.utc).isoformat()}
    return receipt,[{'program':program,**{k:r[k] for k in FIELDS[1:]}} for r in rows]

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=Path(__file__).parent/'archive')
    ap.add_argument('--cache',type=Path,help='optional directory for reusable successful per-program responses')
    args=ap.parse_args()
    path=args.out
    with (path/'visits.csv').open() as f:
        programs=sorted({int(r['program']) for r in csv.DictReader(f)})
    def get_cached(pid):
        if args.cache is None:
            return papers(pid)
        args.cache.mkdir(exist_ok=True, parents=True)
        cached=args.cache/f'{pid}.json'
        if cached.exists():
            return json.loads(cached.read_text())
        result=papers(pid)
        cached.write_text(json.dumps(result))
        return result
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        resp=list(pool.map(get_cached,programs))
    with (path/'program-papers.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS,lineterminator='\n');writer.writeheader()
        for _,rows in resp:
            writer.writerows(rows)
    (path/'papertrack-receipt.json').write_text(json.dumps({'endpoint':ENDPOINT,'program_queries':[r for r,_ in resp]},indent=2)+'\n')
