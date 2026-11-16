"""Sample GOES-19 550 nm aerosol optical depth at overpass and fixed UTC hours at paired PACE coordinates."""
import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path
import sys
from urllib.request import urlopen

from osgeo import gdal, osr

ROOT=Path(__file__).resolve().parent
KEYS=[
    'ABI-L2-AODC/2025/133/17/OR_ABI-L2-AODC-M6_G19_s20251331721177_e20251331723550_c20251331727595.nc',
    'ABI-L2-AODC/2025/133/18/OR_ABI-L2-AODC-M6_G19_s20251331801177_e20251331803550_c20251331807431.nc',
    'ABI-L2-AODC/2025/137/18/OR_ABI-L2-AODC-M6_G19_s20251371801181_e20251371803553_c20251371807375.nc',
    'ABI-L2-AODC/2025/137/18/OR_ABI-L2-AODC-M6_G19_s20251371806181_e20251371808553_c20251371812290.nc',
    'ABI-L2-AODC/2025/141/17/OR_ABI-L2-AODC-M6_G19_s20251411711173_e20251411713546_c20251411716594.nc',
    'ABI-L2-AODC/2025/141/18/OR_ABI-L2-AODC-M6_G19_s20251411801173_e20251411803546_c20251411806546.nc',
    'ABI-L2-AODC/2025/141/18/OR_ABI-L2-AODC-M6_G19_s20251411846173_e20251411848546_c20251411851247.nc',
]
BASE='https://noaa-goes19.s3.amazonaws.com/'
def sample(item):
    key,points=item
    scan=datetime.strptime(key.split("_s")[1][:13],"%Y%j%H%M%S")
    day=scan.date().isoformat()
    path=ROOT/'cache'/Path(key).name
    path.parent.mkdir(exist_ok=True)
    if not path.exists():
        with urlopen(BASE+key,timeout=120) as r: path.write_bytes(r.read())
    gdal.UseExceptions()
    ds=gdal.Open(f'NETCDF:"{path}":AOD')
    flags=gdal.Open(f'NETCDF:"{path}":DQF')
    band=ds.GetRasterBand(1)
    src=osr.SpatialReference()
    src.ImportFromEPSG(4326)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    proj=osr.CoordinateTransformation(src,ds.GetSpatialRef())
    inverse=gdal.InvGeoTransform(ds.GetGeoTransform())
    results=[]
    for p in points:
        if day not in (p['start'],p['end']): continue
        lon,lat=float(p['lon']),float(p['lat'])
        x,y,_=proj.TransformPoint(lon,lat)
        col,row=[int(a) for a in gdal.ApplyGeoTransform(inverse,x,y)]
        assert 0<=col<ds.RasterXSize and 0<=row<ds.RasterYSize
        value=int(band.ReadAsArray(col,row,1,1)[0,0])
        dqf=int(flags.GetRasterBand(1).ReadAsArray(col,row,1,1)[0,0])
        aod=None if dqf==3 or value==band.GetNoDataValue() else round(value*band.GetScale()+band.GetOffset(),6)
        results.append(dict(start=p['start'],end=p['end'],date=day,scan_utc=scan.isoformat()+'Z',lat=lat,lon=lon,col=col,row=row,dqf=dqf,aod_550=aod))
    receipt=dict(date=day,scan_utc=scan.isoformat()+'Z',bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),url=BASE+key)
    return results,receipt

def main():
    with open(sys.argv[1]) as f:
        points=[r for r in csv.DictReader(f) if r['region']=='florida_atlantic']
    results=[sample((key,points)) for key in KEYS]
    rows=[p for r,_ in results for p in r]
    with (ROOT/'goes19_aod.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    (ROOT/'aod_receipts.json').write_text(json.dumps([r for _,r in results],indent=2)+'\n')

if __name__=='__main__': main()
