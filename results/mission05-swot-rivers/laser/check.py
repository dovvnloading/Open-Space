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
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent


def read_observations(record, cache, *, offline=False):
    path = cache / f"{record['date']}_{record['rgt']}.json"
    if not offline and not path.exists():
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
    parser.add_argument('cache', type=Path, nargs='?', help='external cache for the 85-date search')
    parser.add_argument('--event-only', action='store_true',
                        help='verify and print the archived March 2024 comparison offline')
    parser.add_argument('--out', type=Path, default=BASE)
    args = parser.parse_args()
    if args.event_only == (args.cache is not None):
        parser.error('choose --event-only or supply a cache directory')
    with (BASE/'requests.csv').open() as f:
        records = list(csv.DictReader(f))
    if args.event_only:
        records = [r for r in records if (r['date'], r['rgt']) == ('2024-03-13', '1303')]
        if len(records) != 1:
            raise ValueError('Expected one March 13, 2024 laser receipt')
        rows = read_observations(records[0], BASE, offline=True)
    else:
        args.cache.mkdir(parents=True, exist_ok=True)
        with ThreadPoolExecutor(max_workers=6) as pool:
            rows = [p for group in pool.map(lambda r: read_observations(r, args.cache), records)
                    for p in group]
    positions = pd.DataFrame(rows, columns=['date','rgt','beam','lat','lon','z'])
    series = []
    congo = BASE.parent/'congo'
    for receipt in json.loads((congo/'node-manifest.json').read_text()):
        file = congo/receipt['filename']
        raw = file.read_bytes()
        if hashlib.sha256(raw).hexdigest() != receipt['sha256']:
            raise ValueError(f'Changed SWOT node response: {file.name}')
        series.append(pd.read_csv(io.StringIO(json.loads(raw)['results']['csv'])))
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
    matches = pd.DataFrame(results)
    if args.event_only:
        event = matches[matches.radius_km.eq(1.0)].set_index('node_id')
        lower_id, upper_id = 13217000060451, 13217000090231
        expected = pd.read_csv(BASE/'node_matches.csv')
        expected = expected[expected.laser_day.eq('2024-03-13') & expected.radius_km.eq(1.0)].set_index('node_id')
        if set(event.index) != {lower_id, upper_id} or set(expected.index) != set(event.index):
            raise ValueError('Expected both March 2024 node matches at a 1 km radius')
        fields = ['n_segments', 'n_beams', 'laser_median_m', 'swot_m']
        np.testing.assert_allclose(event[fields].sort_index(), expected[fields].sort_index(),
                                   rtol=0, atol=0.00000051)
        lower, upper = event.loc[lower_id], event.loc[upper_id]
        laser_drop = upper.laser_median_m - lower.laser_median_m
        swot_drop = upper.swot_m - lower.swot_m
        print(f'13 March 2024 ATL13 drop: {laser_drop:.6f} m '
              f'({lower.n_segments} downstream and {upper.n_segments} upstream segments)')
        print(f'10 March 2024 SWOT drop: {swot_drop:.6f} m')
        print(f'Laser minus SWOT drop: {laser_drop - swot_drop:.6f} m')
    else:
        matches.to_csv(args.out/'node_matches.csv',index=False,float_format='%.6f')


if __name__ == '__main__':
    main()
