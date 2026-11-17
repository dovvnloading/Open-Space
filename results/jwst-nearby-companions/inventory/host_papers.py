"""Find object-linked imaging/companion papers for the GCNS-matched target inventory."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE=Path(__file__).resolve().parent
URL='https://simbad.cfa.harvard.edu/simbad/sim-tap/sync'
KEYWORDS=" OR ".join(f"title LIKE '{pattern}'" for pattern in ("%maging%","%oronagraph%","%irect %etection%","%igh%ontrast%"))
OBJECTS=" OR ".join(f"title LIKE '{pattern}'" for pattern in ("%ompanion%","%lanet%","%rown %warf%","%ubstellar%","%ub-stellar%"))
def task(job):
    ids,num,refresh=job
    names=','.join("'Gaia DR3 "+x+"'" for x in ids)
    q=("SELECT DISTINCT ident.id, main_id,bibcode,title,doi FROM ident JOIN basic ON ident.oidref=basic.oid "
       "JOIN has_ref ON basic.oid=has_ref.oidref JOIN ref ON ref.oidbib=has_ref.oidbibref "
       "WHERE ident.id IN ("+names+") AND ("+KEYWORDS+") AND ("+OBJECTS+")")
    args=dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC='20000',QUERY=q)
    d=HERE/'archive/host-literature';d.mkdir(exist_ok=True)
    path=d/f'page-{num:02d}.csv.gz';meta=d/f'page-{num:02d}.json'
    if refresh:
        with urlopen(URL+'?'+urlencode(args),timeout=180) as f:raw=f.read()
        if not raw.startswith(b'id,main_id,bibcode,title,doi'):raise RuntimeError(raw[:500])
        path.write_bytes(gzip.compress(raw,mtime=0))
        meta.write_text(json.dumps(dict(endpoint=URL,params=args,sha256=hashlib.sha256(raw).hexdigest(),
                        retrieved_utc=datetime.now(timezone.utc).isoformat()),indent=2)+'\n')
    raw=gzip.decompress(path.read_bytes())
    if hashlib.sha256(raw).hexdigest()!=json.loads(meta.read_text())['sha256']:raise ValueError(path)
    out=list(csv.DictReader(io.StringIO(raw.decode())))
    if len(out)>=20000:raise ValueError("truncation")
    return out

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--refresh',action='store_true')
    a=parser.parse_args()
    with (HERE/'archive/nearby-fields.csv').open() as f:vis=list(csv.DictReader(f))
    fields={}
    for v in vis:
        fields.setdefault(v['gcns_edr3_source_id'],set()).add(v['targprop'])
    ids=sorted(fields)
    with ThreadPoolExecutor(3) as pool:
        matches=[r for subset in pool.map(task,[(ids[i:i+24],i//24,a.refresh) for i in range(0,len(ids),24)]) for r in subset]
    with (HERE/'archive/target-search-status.csv').open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['gcns_edr3_source_id','target_aliases','title_selected_papers','interpretation','title_selected_JWST_papers','bibcodes'])
        for id in ids:
            chosen=[r for r in matches if r['id'].split()[-1]==id]
            jw=[r for r in chosen if 'JWST' in r['title'].upper() or 'JAMES WEBB' in r['title'].upper()]
            w.writerow([id,' | '.join(sorted(fields[id])),len(chosen),
                'object-linked leads; verify target use in paper' if chosen else 'no title-selected reference for this DR3 identifier; prior search status open',
                len(jw),' | '.join(sorted({r['bibcode'] for r in chosen}))])
    print(len(ids),len(matches),len({r['id'] for r in matches}),len({r['bibcode'] for r in matches}))
if __name__=='__main__':main()
