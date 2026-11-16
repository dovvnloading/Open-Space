"""Check independent ICESat-2 ATL13 elevations near six SWOT nodes at Malebo Pool.

Run with an external cache directory. The committed, date/track-specific requests
are a frozen manifest for replay; each JSON payload is verified before use.
"""
import argparse
import csv
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent


def download(record, cache):
    path = cache / f"{record['date']}_{record['rgt']}.json"
    if not path.exists():
        with urlopen(record['atl13_url'],timeout=60) as response:
            path.write_bytes(response.read())
    body = path.read_bytes()
    if hashlib.sha256(body).hexdigest() != record['response_sha256']:
        raise ValueError(f'Changed response: {path.name}')
    response = json.loads(body)
    return [(record['date'],record['rgt'],series['beam'],lat,lon,z)
            for series in response['series'] for lat,lon,z in series['lat_lon_elev']]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache',type=Path)
    parser.add_argument('--out',type=Path,default=BASE)
    args = parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    with (BASE/'requests.csv').open() as f:
        records = list(csv.DictReader(f))
    with ThreadPoolExecutor(max_workers=6) as pool:
        rows = [p for group in pool.map(lambda r:download(r,args.cache),records) for p in group]
    positions = pd.DataFrame(rows,columns=['date','rgt','beam','lat','lon','z'])
    series = []
    for file in sorted((BASE.parent/'congo/raw').glob('node-*.json')):
        series.append(pd.read_csv(io.StringIO(json.loads(file.read_text())['results']['csv'])))
    if len(series)!=6:
        raise ValueError('Expected six archived SWOT node series')
    nodes = pd.concat(series)
    nodes = nodes[(nodes.node_q<=1)&nodes.wse_u.between(0,.30,inclusive='right')&(nodes.wse>0)].copy()
    nodes['date'] = pd.to_datetime(nodes.time_str,utc=True).dt.tz_localize(None).dt.normalize()
    results = []
    for day,v in positions.groupby('date'):
        deltas = (nodes.date-pd.Timestamp(day)).dt.days.abs()
        for node, group in nodes[deltas<=7].groupby('node_id'):
            # Acquisition chosen before inspecting or averaging laser heights.
            sat = group.iloc[np.argmin((group.date-pd.Timestamp(day)).abs().to_numpy())]
            distance = np.hypot((v.lon-sat.lon)*111.32*np.cos(np.deg2rad(sat.lat)),(v.lat-sat.lat)*110.57)
            for radius in (.5,1.0):
                pts = v.loc[distance<=radius]
                if len(pts)<3:
                    continue
                results.append(dict(laser_day=day,swot_day=sat.date.date(),rgt=pts.rgt.iloc[0],
                                    node_id=node,pass_id=sat.pass_id,cycle_id=sat.cycle_id,
                                    radius_km=radius,n_segments=len(pts),n_beams=pts.beam.nunique(),
                                    min_km=distance.min(),median_km=distance.loc[pts.index].median(),
                                    laser_median_m=pts.z.median(),laser_min_m=pts.z.min(),
                                    laser_max_m=pts.z.max(),swot_m=sat.wse,swot_u_m=sat.wse_u,
                                    laser_minus_swot_m=pts.z.median()-sat.wse))
    pd.DataFrame(results).to_csv(args.out/'node_matches.csv',index=False,float_format='%.6f')


if __name__ == '__main__':
    main()
