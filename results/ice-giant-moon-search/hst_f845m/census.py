"""Fetch the compact WFC3 F845M sample for reuse across HST screening passes."""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
from pathlib import Path
import sys
import time
from urllib.request import urlopen
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'hst16995'))
import reduce


def download(row,cache):
    n=row['rootname'];fn=n+'_flt.fits';url=f'https://stpubdata.s3.amazonaws.com/hst/public/{n[:4]}/{n}/{fn}'
    path=cache/fn
    if path.exists(): return
    for trial in range(4):
        try:
            with urlopen(url,timeout=90) as response:raw=response.read()
            if not raw.startswith(b'SIMPLE  =') or len(raw)%2880:raise ValueError(fn)
            path.write_bytes(raw);return
        except Exception:
            if trial==3:raise
            time.sleep(2**trial)


def acquire(cache,freeze):
    inventory=list(csv.DictReader((ROOT.parent/'manifests/observations.csv').open()))
    rows=[]
    for r in sorted(inventory,key=lambda d:float(d['t_min'])):
        if r['obs_collection']!='HST' or r['instrument_name']!='WFC3/UVIS' or r['filters']!='F845M':continue
        n=r['obs_id'];rows.append(dict(target=r['planet'].upper(),obsid=r['obsid'],rootname=n,visit=n[4:6],mjd_start=r['t_min'],mjd_end=r['t_max'],exposure_s=r['t_exptime'],uri='mast:HST/product/'+n+'_flt.fits',proposal=r['proposal_id']))
    if freeze:reduce.write_csv(ROOT/'selection.csv',rows)
    else:rows=list(csv.DictReader((ROOT/'selection.csv').open()))
    cache.mkdir(exist_ok=True,parents=True)
    with ThreadPoolExecutor(4) as ex:list(ex.map(lambda r:download(r,cache),rows))
    hashes=[]
    for r in rows:
        name=r['rootname']+'_flt.fits'
        path=cache/name
        hashes.append(dict(file=name,bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    if freeze:
        reduce.write_csv(ROOT/'inputs.csv',hashes)
    else:
        frozen=list(csv.DictReader((ROOT/'inputs.csv').open()))
        if frozen != [{k:str(v) for k,v in r.items()} for r in hashes]:
            raise ValueError('F845M input checksum or size changed')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('cache',type=Path);ap.add_argument('--freeze',action='store_true')
    a=ap.parse_args();acquire(a.cache,a.freeze)
