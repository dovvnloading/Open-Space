"""Fetch the 11x11-pixel TESS observations for every covered sector in the inventory."""
import concurrent.futures
import json
import sys
import time
import zipfile
from pathlib import Path
from io import BytesIO
from urllib.parse import urlencode
from urllib.request import urlopen

HERE=Path(__file__).resolve().parent
BASE='https://mast.stsci.edu/tesscut/api/v0.1/'
RA,DEC=21.06275,-73.32594444

def url(s):
    return BASE+'astrocut?'+urlencode(dict(ra=RA,dec=DEC,sector=s,x=11,y=11,units='px'))

def fetch(s):
    path=HERE/'raw'/f's{s:03}.fits'
    if path.exists(): return
    for attempt in range(3):
        try:
            with urlopen(url(s),timeout=140) as r: blob=r.read()
            with zipfile.ZipFile(BytesIO(blob)) as z:
                entries=[n for n in z.namelist() if n.endswith('.fits')]
                if len(entries)!=1: raise ValueError(entries)
                raw=z.read(entries[0])
            temp=path.with_suffix('.part');temp.write_bytes(raw);temp.rename(path)
            print(s,len(raw),flush=True);return
        except Exception as e:
            print(s,type(e).__name__,str(e)[:90],flush=True)
            if attempt==2: raise
            time.sleep(4)

if __name__=='__main__':
    (HERE/'raw').mkdir(exist_ok=True)
    inventory=HERE/'sectors.json'
    if not inventory.exists():
        with urlopen(BASE+'sector?'+urlencode(dict(ra=RA,dec=DEC)),timeout=60) as r:inventory.write_bytes(r.read())
    sectors=[int(row['sector']) for row in json.loads(inventory.read_text())['results']]
    with concurrent.futures.ThreadPoolExecutor(3) as workers:
        list(workers.map(fetch,sectors))
