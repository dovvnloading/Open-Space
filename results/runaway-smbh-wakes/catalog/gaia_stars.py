"""Download bright Gaia DR3 stars around the four Q1 field groups through the Heidelberg Gaia TAP service."""
import csv
import gzip
import hashlib
import io
import os
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE=Path(__file__).resolve().parent
CIRCLES={'EDF-N':(269.8,66,6),'EDF-S':(61.1,-48.5,6),'EDF-F':(53,-28,4.5),'LDN1641':(85.8,-8.5,2)}

def run(item):
    label,(ra,dec,r)=item
    query=(f'SELECT source_id,ra,dec,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE '
           f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{r})) AND phot_g_mean_mag<16")
    url='https://gaia.ari.uni-heidelberg.de/tap/sync?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC=1000000,QUERY=query))
    with urlopen(url,timeout=240) as response:raw=response.read()
    assert raw.startswith(b'source_id,'),raw[:250]
    rows=list(csv.DictReader(io.StringIO(raw.decode())))
    assert len(rows)<1000000
    cache=Path(os.environ.get('WAKE_GAIA_CACHE',HERE/'gaia_cache'));cache.mkdir(exist_ok=True,parents=True)
    file=cache/('gaia_'+label+'.csv.gz')
    with file.open('wb') as f,gzip.GzipFile(fileobj=f,mode='wb',mtime=0) as g:g.write(raw)
    return dict(field=label,url=url,rows=len(rows),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())

if __name__=='__main__':
    with ThreadPoolExecutor(4) as p: out=list(p.map(run,CIRCLES.items()))
    (HERE/'gaia_queries.json').write_text(json.dumps(out,indent=2)+'\n')
