"""Measure Sentinel-2 classified open water in fixed SWORD reach corridors.

Run with Python 3, NumPy, GDAL (osgeo) and network access.
The cache stores the original SCL rasters and the catalog responses outside Git.
"""
import argparse
import calendar
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
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
BASE = Path(__file__).resolve().parent
SWORD = BASE.parent / 'zambezi/data/sword_v16_selected.geojson'
SEARCH = 'https://earth-search.aws.element84.com/v1/search'
# Metres either side of the fixed SWORD v16 centreline. The half-widths
# exceed half the widest SWORD node width in each study reach.
SITES = {'12297000161': 800, '12279000221': 1400}


def cached(path, url):
    if not path.exists():
        for attempt in range(3):
            try:
                with urlopen(url, timeout=50) as r:
                    data = r.read()
                path.write_bytes(data)
                break
            except OSError:
                if attempt == 2:
                    raise
                time.sleep(2 ** attempt)
    return path.read_bytes()


def search(task):
    site, month, cache = task
    props = site['properties']
    x, y = props['x'], props['y']
    end = month.replace(day=calendar.monthrange(month.year, month.month)[1])
    url = SEARCH + '?' + urlencode(dict(collections='sentinel-2-l2a', bbox=f'{x-.001},{y-.001},{x+.001},{y+.001}',
                                        limit=100, datetime=f'{month}T00:00:00Z/{end}T23:59:59Z'))
    body = cached(cache / f"{props['reach_id']}_{month}.json", url)
    data = json.loads(body)
    assert len(data['features']) == data['numberMatched']
    # Up to two dates per month. Use tile cloud metadata, prior to inspecting local water.
    by_day = {}
    for item in sorted(data['features'], key=lambda i: (i['properties']['eo:cloud_cover'], i['id'])):
        by_day.setdefault(item['properties']['datetime'][:10], item)
    result = []
    for item in list(by_day.values())[:2]:
        result.append((site, item, cache, url, hashlib.sha256(body).hexdigest()))
    return result


def measure(task):
    site, item, cache, url, catalog_hash = task
    reach = site['properties']['reach_id']
    scene = item['id']
    href = item['assets']['scl']['href']
    path = cache / (scene + '.tif')
    raw = cached(path, href)
    src = gdal.Open(str(path))
    srs = src.GetSpatialRef()
    wgs = osr.SpatialReference()
    wgs.ImportFromEPSG(4326)
    wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    line = ogr.CreateGeometryFromJson(json.dumps(site['geometry']))
    line.Transform(osr.CoordinateTransformation(wgs, srs))
    corridor = line.Buffer(SITES[reach])
    xmin, xmax, ymin, ymax = corridor.GetEnvelope()
    gt = src.GetGeoTransform()
    assert gt[1] == 20 and gt[5] == -20 and gt[2] == gt[4] == 0
    left, top = int(np.floor((xmin-gt[0])/20)), int(np.floor((gt[3]-ymax)/20))
    right, bottom = int(np.ceil((xmax-gt[0])/20)), int(np.ceil((gt[3]-ymin)/20))
    # An overlapping tile that does not contain the whole corridor cannot be compared.
    rec = dict(reach_id=reach, date=item['properties']['datetime'][:10], scene=scene,
               baseline=item['properties'].get('s2:processing_baseline',''),
               tile_cloud_pct=item['properties']['eo:cloud_cover'], query_url=url,
               query_sha256=catalog_hash, scl_url=href, scl_sha256=hashlib.sha256(raw).hexdigest(),
               buffer_m=SITES[reach], pixels='', clear_pixels='', water_pixels='', clear_fraction='',
               open_water_km2='', max_open_water_km2='', accepted=0, reason='partial_tile')
    if left<0 or top<0 or right>src.RasterXSize or bottom>src.RasterYSize:
        return rec
    nx, ny = right-left, bottom-top
    arr = src.GetRasterBand(1).ReadAsArray(left,top,nx,ny)
    mask_ds = gdal.GetDriverByName('MEM').Create('',nx,ny,1,gdal.GDT_Byte)
    mask_ds.SetProjection(src.GetProjection())
    mask_ds.SetGeoTransform((gt[0]+left*20,20,0,gt[3]-top*20,0,-20))
    mem = ogr.GetDriverByName('Memory').CreateDataSource('')
    layer = mem.CreateLayer('corridor',srs,ogr.wkbPolygon)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(corridor)
    layer.CreateFeature(feature)
    gdal.RasterizeLayer(mask_ds,[1],layer,burn_values=[1])
    vals = arr[mask_ds.ReadAsArray().astype(bool)]
    clear = np.isin(vals,[4,5,6]).sum()
    water = (vals==6).sum()
    accepted = clear / vals.size >= .98
    rec.update(pixels=vals.size,clear_pixels=clear,water_pixels=water,
               clear_fraction=round(clear/vals.size,6),open_water_km2=round(water*.0004,5),
               max_open_water_km2=round((water+vals.size-clear)*.0004,5),
               accepted=int(accepted),reason='accepted' if accepted else 'obscured')
    return rec


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cache',type=Path)
    p.add_argument('--out',type=Path,default=BASE/'scenes.csv')
    args = p.parse_args()
    args.cache.mkdir(exist_ok=True,parents=True)
    sites = [f for f in json.loads(SWORD.read_text())['features'] if f['properties']['reach_id'] in SITES]
    month, dates = date(2023,7,1), []
    while month <= date(2025,4,1):
        dates.append(month)
        month = (month+timedelta(days=32)).replace(day=1)
    with ThreadPoolExecutor(max_workers=8) as pool:
        selection = list(pool.map(search, ((s,m,args.cache) for s in sites for m in dates)))
        rows = list(pool.map(measure,[t for group in selection for t in group]))
    with args.out.open('w') as out:
        writer = csv.DictWriter(out,fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f'{len(rows)} scenes; {sum(r["accepted"] for r in rows)} pass the local clear-area test')


if __name__ == '__main__':
    main()
