"""Localize the January–July 2024 stage response along all SWORD nodes of two reaches.

Query the two pass-111 windows separately. Build a fixed-node spatial profile and
evaluate a second track on separate dates. Raw response hashes and URLs are recorded.
"""
import argparse
import csv
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import time
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
# One query per node retains the four windows via a single acquisition request.
WINDOWS = {'target_high':'2024-01-08','target_low':'2024-07-13',
           'check_high':'2023-12-13','check_low':'2024-08-20'}
NODES = [f'{rid[:-1]}{idx:03}1' for rid,n in [('13217000061',60),('13217000091',46)] for idx in range(1,n+1)]
URL = 'https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries'
FIELDS = 'node_id,reach_id,time_str,cycle_id,pass_id,lat,lon,wse,wse_u,node_q,node_q_b,crid,granuleUR'


def source(task):
    node,cache = task
    path = cache/f'{node}.json'
    url = URL+'?'+urlencode(dict(feature='Node',feature_id=node,output='csv',start_time='2023-12-12T00:00:00Z',
                    end_time='2024-08-21T00:00:00Z',collection_name='SWOT_L2_HR_RiverSP_2.0',fields=FIELDS))
    if not path.exists():
        for attempt in range(3):
            try:
                with urlopen(url,timeout=50) as f:
                    raw=f.read()
                parsed=json.loads(raw)
                if parsed['status']!='200 OK':
                    raise RuntimeError(parsed['status'])
                path.write_bytes(raw)
                break
            except (OSError,ValueError,RuntimeError):
                if attempt==2:
                    raise
                time.sleep(1+attempt)
    raw=path.read_bytes()
    payload=json.loads(raw)['results']['csv']
    rows=list(csv.DictReader(io.StringIO(payload)))
    selected=[r for r in rows if r['time_str'][:10] in WINDOWS.values() or
              ((r['cycle_id'],r['pass_id']) in {('009','111'),('018','111'),('007','570'),('019','570')})]
    return dict(node_id=node,url=url,rows=len(rows),response_sha256=hashlib.sha256(raw).hexdigest(),
                csv_sha256=hashlib.sha256(payload.encode()).hexdigest()),selected


def acquire(cache):
    cache.mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(max_workers=6) as pool:
        data=list(pool.map(source,[(node,cache) for node in NODES]))
    pd.DataFrame(x[0] for x in data).to_csv(BASE/'manifest.csv',index=False)
    pd.DataFrame(r for _,rows in data for r in rows).to_csv(BASE/'measurements.csv',index=False)


def analyze():
    a=pd.read_csv(BASE/'measurements.csv').replace(-999999999999.,np.nan)
    a['day']=a.time_str.str[:10]
    tables=[]
    for pair in ('target','check'):
        g=a[a.day.eq(WINDOWS[pair+'_high'])].merge(a[a.day.eq(WINDOWS[pair+'_low'])],on='node_id',
                                                  suffixes=('_high','_low'),validate='one_to_one')
        for phase in ('high','low'):
            g=g[g['node_q_'+phase].isin([0,1])&g['wse_u_'+phase].between(0,.30,inclusive='right')]
        g=g.dropna(subset=['wse_high','wse_low'])
        g['rise_m']=g.wse_high-g.wse_low
        g['run']=pair
        g['node_idx']=g.node_id.astype(str).str[-4:-1].astype(int)
        g['shift_m']=np.hypot((g.lat_high-g.lat_low)*110574,
                    (g.lon_high-g.lon_low)*111320*np.cos(np.radians(g.lat_low)))
        tables.append(g[['run','reach_id_high','node_id','node_idx','rise_m','lat_low','lon_low','shift_m']])
    profile=pd.concat(tables)
    profile.to_csv(BASE/'profile.csv',index=False,float_format='%.8f')
    # Statistics for all non-overlapping 5-node bins, ordered downstream to upstream.
    profile['bin_start']=(profile.node_idx-1)//5*5+1
    bins=profile.groupby(['run','reach_id_high','bin_start']).agg(n=('rise_m','count'),rise_m=('rise_m','median'),
                    lat=('lat_low','median'),lon=('lon_low','median'),median_shift_m=('shift_m','median')).reset_index()
    bins.to_csv(BASE/'bins.csv',index=False,float_format='%.8f')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--retrieve',type=Path,metavar='CACHE')
    args=p.parse_args()
    if args.retrieve:
        acquire(args.retrieve)
    analyze()
