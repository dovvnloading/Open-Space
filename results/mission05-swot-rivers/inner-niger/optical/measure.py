"""Measure Sentinel-2 classified open water at the delta outlet and Lake Débo outlet.

Run with system Python (GDAL and NumPy). Store complete source SCL rasters outside Git.
Each target is a fixed 4 x 4 km UTM square centered on the station vicinity.
"""
import argparse
import csv
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import date, timedelta
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()
ROOT = Path(__file__).resolve().parent
SITES = {'dire': (-3.3563, 16.2666, '30QVD'), 'akka': (-4.23, 15.39, '30PUC')}
SEARCH = 'https://earth-search.aws.element84.com/v1/search'


def acquire(url, path):
    if not path.exists():
        for attempt in range(3):
            try:
                with urlopen(url, timeout=50) as response:
                    body = response.read()
                path.write_bytes(body)
                break
            except OSError:
                if attempt == 2:
                    raise
                time.sleep(2 ** attempt)
    return path.read_bytes()


def measure(task):
    station, season, anchor, cache = task
    lon, lat, tile = SITES[station]
    args = dict(collections='sentinel-2-l2a', bbox=f'{lon-.001},{lat-.001},{lon+.001},{lat+.001}',
                datetime=f'{anchor-timedelta(days=14)}T00:00:00Z/{anchor+timedelta(days=14)}T23:59:59Z',
                limit=100)
    query = SEARCH + '?' + urlencode(args)
    raw = acquire(query, cache / f'{station}-{anchor}.json')
    catalog = json.loads(raw)
    if len(catalog['features']) != catalog['numberMatched']:
        raise ValueError('Search requires pagination')
    # One fixed tile, platform A/B. Give the nearest clear-tile date precedence.
    scenes = [s for s in catalog['features'] if s['id'].split('_')[1] == tile and
              s['id'].startswith(('S2A_', 'S2B_')) and s['properties']['eo:cloud_cover'] <= 15]
    scenes.sort(key=lambda s: (abs((date.fromisoformat(s['properties']['datetime'][:10])-anchor).days),
                               s['properties']['eo:cloud_cover'], s['id']))
    rec = dict(station=station, season=season, anchor=str(anchor), query=query,
               catalog_sha256=hashlib.sha256(raw).hexdigest(), candidates=len(scenes), status='no_scene')
    if not scenes:
        return rec
    item = scenes[0]
    asset = item['assets']['scl']['href']
    path = cache / (item['id']+'.tif')
    raw = acquire(asset, path)
    ds = gdal.Open(str(path))
    geo = osr.SpatialReference()
    geo.ImportFromEPSG(4326)
    geo.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    projected = ds.GetSpatialRef()
    x, y, _ = osr.CoordinateTransformation(geo, projected).TransformPoint(lon, lat)
    x0, dx, _, y0, _, dy = ds.GetGeoTransform()
    assert dx == 20 and dy == -20 and projected.GetAuthorityCode(None) == '32630'
    col, row = int((x-x0)//20)-100, int((y0-y)//20)-100
    assert 0 <= col <= ds.RasterXSize-200 and 0 <= row <= ds.RasterYSize-200
    vals = ds.GetRasterBand(1).ReadAsArray(col,row,200,200)
    good = np.isin(vals,[4,5,6]); water = vals==6
    rec.update(status='usable' if good.mean()>=.98 else 'obscured', scene=item['id'],
               utc=item['properties']['datetime'], processing_baseline=item['properties'].get('s2:processing_baseline'),
               asset=asset, raster_bytes=len(raw), raster_sha256=hashlib.sha256(raw).hexdigest(),
               pixel_window_sha256=hashlib.sha256(vals.tobytes()).hexdigest(),
               x_left=x0+20*col,y_top=y0-20*row,clear_pixels=int(good.sum()),water_pixels=int(water.sum()),
               area_lower_km2=round(float(water.sum())*.0004,6),
               area_upper_km2=round(float((water | ~good).sum())*.0004,6))
    return rec


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cache', type=Path)
    a=p.parse_args(); a.cache.mkdir(exist_ok=True, parents=True)
    tasks=[(site, yr, day, a.cache) for site in SITES for yr in range(2015,2026)
           for day in (date(yr,12,15),date(yr+1,1,15))]
    # Four matched hydraulic visits; distinct from the annual crest windows.
    event_days = [date(2024,7,22),date(2025,1,26),date(2025,7,12),date(2026,1,16)]
    tasks.extend(('dire', 'matched_visit', day, a.cache) for day in event_days)
    tasks = list(dict.fromkeys(tasks))
    with ThreadPoolExecutor(max_workers=5) as pool:
        records=list(pool.map(measure,tasks))
    fields=list(dict.fromkeys(key for rec in records for key in rec))
    with (ROOT/'windows.csv').open('w') as f:
        writer=csv.DictWriter(f,fields,lineterminator='\n');writer.writeheader();writer.writerows(records)
    print(len(records),'windows;',sum(r['status']=='usable' for r in records),'usable')


if __name__=='__main__':
    main()
