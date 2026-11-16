"""Download and classify 20 m Sentinel-2 L2A SCL masks around Mopti's primary SWOT reach.

Requires Python 3, GDAL's Python bindings, and NumPy. Run in the optical directory.
Complete (source) files are checked by digest, then a lossless 500x500 window is kept.
"""
import csv
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from hashlib import sha256
import json
from pathlib import Path
import time
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from osgeo import gdal

BASE = Path(__file__).resolve().parent
CACHE = BASE.parent.parent / 'cache'
BBOX = '-4.202,14.532,-4.201,14.533'
# 10 x 10 km grid around the reach midpoint, rounded to the closest 20 m pixel edge.
WEST, NORTH, SIZE = 365560, 1611980, 500
SEARCH = 'https://earth-search.aws.element84.com/v1/search?'

def get(url):
    for attempt in range(3):
        try:
            with urlopen(url, timeout=50) as response:
                return response.read()
        except OSError:
            if attempt == 2:
                raise
            time.sleep(0.5 + attempt)

def write_table(name, rows):
    with (BASE / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def observe_year(year):
    target = date(year,11,8)
    params = {'collections':'sentinel-2-l2a', 'bbox': BBOX,
              'datetime':f'{year}-11-01T00:00:00Z/{year}-11-15T23:59:59Z','limit':100}
    query = SEARCH + urlencode(params)
    catalog_path = CACHE / f'mopti_stac_{year}.json'
    if not catalog_path.exists():
        catalog_path.write_bytes(get(query))
    catalog = json.loads(catalog_path.read_bytes())
    assert catalog['numberMatched'] == len(catalog['features']), 'Incomplete candidate search'
    scenes = [s for s in catalog['features'] if s['properties']['grid:code']=='MGRS-30PUB'
              and '_R008_' in s['properties']['s2:product_uri']]
    def rank(s):
        return abs((date.fromisoformat(s['properties']['datetime'][:10])-target).days),s['properties']['datetime'][:10],-float(s['properties']['s2:processing_baseline']),s['id']
    scenes.sort(key=rank)
    attempts = []
    for s in scenes:
        asset=s['assets']['scl']['href']
        raw = CACHE / (s['id'] + '_SCL.tif')
        if not raw.exists():
            raw.write_bytes(get(asset))
        ds = gdal.Open(str(raw))
        x, scale, _, y, _, scale_y = ds.GetGeoTransform()
        assert scale == 20 and scale_y == -20 and ds.GetSpatialRef().GetAuthorityCode(None)=='32630'
        xoff, yoff = (WEST-x)/scale, (NORTH-y)/scale_y
        assert int(xoff)==xoff and int(yoff)==yoff
        a = ds.ReadAsArray(int(xoff),int(yoff),SIZE,SIZE)
        clear = np.isin(a, [4,5,6])
        water = a == 6
        p = s['properties']
        row = {'year':year,'scene':s['id'],'time_utc':p['datetime'],'relative_orbit':8,
               'processing_baseline':p.get('s2:processing_baseline'),'source_url':asset,
               'source_sha256':sha256(raw.read_bytes()).hexdigest(),'water_pixels':int(water.sum()),
               'certain_land_or_water_pixels':int(clear.sum()), 'total_pixels':a.size,
               'cloud_shadow_snow_nodata_pixels':int(np.isin(a, [0,1,3,8,9,10,11]).sum()),
               'water_km2_lower_bound':round(water.sum()*.0004,4),
               'water_km2_upper_bound':round((water.sum()+a.size-clear.sum())*.0004,4),
               'stac_search_url':query}
        attempts.append(row.copy())
        if np.mean(np.isin(a, [0,1,3,8,9,10,11])) <= .05:
            chip_path = BASE / 'chips' / (s['id'] + '_mopti_SCL.tif')
            gdal.Translate(str(chip_path),ds,srcWin=[int(xoff),int(yoff),SIZE,SIZE],
                           creationOptions=['COMPRESS=DEFLATE'])
            row['chip_sha256'] = sha256(chip_path.read_bytes()).hexdigest()
            return attempts, row
    raise ValueError(f'No qualifying scene for {year}; candidates {len(scenes)}')

def main():
    gdal.UseExceptions()
    CACHE.mkdir(parents=True,exist_ok=True)
    (BASE/'chips').mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(4) as pool:
        results = list(pool.map(observe_year, range(2018,2026)))
    write_table('candidates_examined.csv', [x for attempts,_ in results for x in attempts])
    write_table('mopti_early_november_water.csv', [row for _,row in results])

if __name__ == '__main__':
    main()
