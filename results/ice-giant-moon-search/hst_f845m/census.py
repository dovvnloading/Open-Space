"""Freeze the complete inventory-selected HST F845M sample and run each short epoch."""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time
from urllib.request import urlopen
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'hst16995'))
import fetch
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
    old=fetch.ROOT;fetch.ROOT=ROOT
    try:
        for trial in range(6):
            try:
                fetch.run(cache,freeze);break
            except Exception:
                if trial==5:raise
                time.sleep(2**trial)
    finally:fetch.ROOT=old


def screen(cache):
    original=reduce.read_csv
    load=lambda n:list(csv.DictReader((ROOT/n).open()))
    eph=load('ephemerides.csv'); hashes=load('inputs.csv');selection=load('selection.csv')
    summary=[]
    for planet in reduce.DOMAINS:
        proposals=sorted({r['proposal'] for r in selection if r['target']==planet.upper()})
        for pid in proposals:
            rows=sorted([r for r in selection if r['target']==planet.upper() and r['proposal']==pid],key=lambda r:float(r['mjd_start']))
            chunks=np.split(np.arange(len(rows)),np.flatnonzero(np.diff([float(r['mjd_start']) for r in rows])>3)+1)
            for indices in chunks:
                part=[rows[k] for k in indices]; names={r['rootname'] for r in part};begin=part[0]['mjd_start'].split('.')[0]
                tag=f'{planet.lower()}-{pid}-{begin}'; visits={r['visit'] for r in part}
                result=dict(epoch=tag,archive_exposures=len(part),archive_visits=len(visits),used_exposures=0,used_visits=0,
                            bank=0,supported=0,screen_hits=0,dropout_hits=0,highest_dropout_rank='',status='under_four_visits')
                if len(visits)>=4:
                    mapping={'selection.csv':part,'inputs.csv':hashes,'ephemerides.csv':[r for r in eph if r['rootname'] in names]}
                    reduce.read_csv=lambda n:mapping[n]
                    try:
                        out=ROOT/tag
                        reduce.run(argparse.Namespace(out=out,cache=cache,planet=planet,no_plants=True))
                        s=json.loads((out/'summary.json').read_text())
                        result.update(used_exposures=s['exposures'],used_visits=s['visits'],bank=s['grid_tracks'],supported=s['supported_tracks'],
                                      screen_hits=s['threshold_crossings'],dropout_hits=s['dropout_crossings'],highest_dropout_rank=round(s['top_dropout_rank'],4),status='screened')
                    finally:reduce.read_csv=original
                summary.append(result);reduce.write_csv(ROOT/'coverage.csv',summary)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('action',choices=['acquire','screen'])
    ap.add_argument('cache',type=Path);ap.add_argument('--freeze',action='store_true');a=ap.parse_args()
    if a.action=='acquire':acquire(a.cache,a.freeze)
    else:screen(a.cache)
