"""Recompute area-weighted July–September CHIRPS v2 precipitation, using GDAL + NumPy."""
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from pathlib import Path
from urllib.request import urlopen
import csv
import gzip
import time

import numpy as np
from osgeo import gdal, ogr, osr

HERE = Path(__file__).resolve().parent
CACHE = HERE.parent.parent / 'cache' / 'chirps'
BASE = 'https://data.chc.ucsb.edu/products/CHIRPS/v2.0/africa_monthly/tifs/'
# This window lies on the CHIRPS 0.05-degree grid and encloses all selected polygons.
WEST, NORTH, COLS, ROWS = -13, 17, 240, 200

def masks():
    vector = ogr.Open(str(HERE / 'mopti_upstream_level6.geojson'))
    ly = vector.GetLayer()
    xmin,xmax,ymin,ymax = ly.GetExtent()
    assert WEST<=xmin and xmax<=WEST+.05*COLS and NORTH-.05*ROWS<=ymin and ymax<=NORTH
    ref=osr.SpatialReference()
    ref.ImportFromEPSG(4326)
    layers=[]
    for upstream_only in (False,True):
        ds=gdal.GetDriverByName('MEM').Create('',COLS,ROWS,1,gdal.GDT_Byte)
        ds.SetGeoTransform((WEST,.05,0,NORTH,0,-.05))
        ds.SetProjection(ref.ExportToWkt())
        ly.SetAttributeFilter('HYBAS_ID <> 1060603890' if upstream_only else None)
        ly.ResetReading()
        gdal.RasterizeLayer(ds, [1], ly, burn_values=[1])
        layers.append(ds.ReadAsArray().astype(bool))
    lat=np.arange(ROWS)*(-.05)+NORTH-.025
    weights=np.broadcast_to(np.cos(np.radians(lat))[:,None],(ROWS,COLS))
    return [(m,weights[m]) for m in layers]

def summarize(param):
    year, month, mask_weights = param
    name=f'chirps-v2.0.{year}.{month:02}.tif.gz'
    p=CACHE/name
    uri=BASE+name
    for attempt in range(3):
        try:
            if not p.exists():
                with urlopen(uri,timeout=45) as r: data=r.read()
                assert data[:2] == bytes.fromhex('1f8b')
                p.write_bytes(data)
            break
        except OSError:
            if attempt==2: raise
            time.sleep(1+attempt)
    compressed=p.read_bytes()
    virtual='/vsimem/'+name[:-3]
    gdal.FileFromMemBuffer(virtual,gzip.decompress(compressed))
    try:
        d=gdal.Open(virtual)
        x,sx,_,y,_,sy=d.GetGeoTransform()
        assert abs(sx-.05)<1e-7 and abs(sy+.05)<1e-7
        a=d.ReadAsArray(round((WEST-x)/sx),round((NORTH-y)/sy),COLS,ROWS)
        row={'year':year,'month':month}
        for (mask,w),key in zip(mask_weights,('rain_mm','upstream_only_mm')):
            vals=a[mask]
            assert np.isfinite(vals).all() and vals.min()>=0
            row[key]=round(float(np.average(vals,weights=w)),5)
        row.update({'source_sha256':sha256(compressed).hexdigest(),'source_bytes':len(compressed),'url':uri})
        return row
    finally:
        gdal.Unlink(virtual)

def table(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)

def main():
    gdal.UseExceptions()
    CACHE.mkdir(parents=True,exist_ok=True)
    mask_weights=masks()
    tasks=[(yr,mo,mask_weights) for yr in range(1991,2026) for mo in (7,8,9)]
    with ThreadPoolExecutor(8) as pool:
        months=list(pool.map(summarize,tasks))
    table(HERE/'chirps_monthly.csv',months)
    years=[]
    for yr in range(1991,2026):
        rows=[r for r in months if r['year']==yr]
        years.append({'flood_year':yr,'july_to_september_mm':round(sum(r['rain_mm'] for r in rows),4),'upstream_only_mm':round(sum(r['upstream_only_mm'] for r in rows),4)})
    table(HERE/'chirps_seasonal.csv',years)

if __name__=='__main__':
    main()
