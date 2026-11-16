"""Extract physical fields at the Florida paired samples, including prior-May controls."""
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path
import re
import sys
import time
from urllib.parse import quote
from urllib.request import urlopen

import numpy as np

ROOT = Path(__file__).resolve().parent
PSL = 'https://psl.noaa.gov/thredds/dodsC/Datasets/'
GEO = 'https://cwcgom.aoml.noaa.gov/thredds/dodsC/OCEAN_GEOSTROPHIC_CURRENTS/CURRENTS.nc'
YEARS = range(2017, 2026)  # First full May of the physical current product is 2017.
FIELDS = ('anom','u_current','v_current','uwnd','vwnd')

def grid(year, field):
    start = date(year,5,1)
    if field == 'anom':
        base = PSL+f'noaa.oisst.v2.highres/sst.day.anom.{year}.nc'
        ref = date(year,1,1)
        y,x=(467,482),(1118,1126)
        epoch,scale=date(1800,1,1),1
    elif field.endswith('current'):
        base = GEO
        ref = date(2016,7,1)
        y,x=(457,474),(498,505)
        epoch,scale=date(1970,1,1),86400
    else:
        base = PSL+f'ncep.reanalysis.dailyavgs/surface_gauss/{field}.10m.gauss.{year}.nc'
        ref = date(year,1,1)
        y,x=(30,33),(148,151)
        epoch,scale=date(1800,1,1),24
    t=(start-ref).days
    query=f'{field}[{t}:1:{t+30}][{y[0]}:1:{y[1]}][{x[0]}:1:{x[1]}]'
    url=base+'.ascii?'+quote(query)
    cache=ROOT/'cache'
    cache.mkdir(exist_ok=True)
    path=cache/f'{year}-may-{field}.txt'
    if path.exists(): raw=path.read_bytes()
    else:
        for attempt in range(4):
            try:
                with urlopen(url,timeout=90) as r: raw=r.read()
                break
            except OSError:
                if attempt==3: raise
                time.sleep(2**attempt)
        path.write_bytes(raw)
    payload=raw.decode().split('---------------------------------------------\n')[1]
    def coords(key):
        block=re.search(rf'{field}\.{key}\[\d+\]\n([^\n]+)',payload)
        return np.fromstring(block.group(1),sep=',')
    lat,lon=coords('lat'),(coords('lon')+180)%360-180
    times=coords('time')/scale
    assert np.array_equal(times,(start-epoch).days+np.arange(31))
    data_rows=re.findall(r'^\[\d+\]\[\d+\], ([^\n]+)',payload,re.M)
    values=np.asarray([np.fromstring(row,sep=',') for row in data_rows]).reshape(31,len(lat),len(lon))
    values[np.abs(values)>1e20]=np.nan
    return (year,field,lat,lon,values,dict(year=year,field=field,bytes=len(raw),
            sha256=hashlib.sha256(raw).hexdigest(),url=url))

def main():
    with open(sys.argv[1]) as fp:
        optical=[p for p in csv.DictReader(fp) if p['region']=='florida_atlantic']
    output=[]; baseline=[];receipts=[]
    with ThreadPoolExecutor(5) as pool:
        results=pool.map(lambda args:grid(*args),[(y,f) for y in YEARS for f in FIELDS])
        for year,field,lat,lon,values,receipt in results:
            receipts.append(receipt)
            for start in ('2025-05-13','2025-05-17'):
                pairs=[r for r in optical if r['start']==start]
                y=np.array([int(np.argmin(np.abs(lat-float(p['lat'])))) for p in pairs])
                x=np.array([int(np.argmin(np.abs(lon-float(p['lon'])))) for p in pairs])
                assert all(abs(lat[i]-float(p['lat']))<=1 for p,i in zip(pairs,y))
                assert all(abs(lon[i]-float(p['lon']))<=1 for p,i in zip(pairs,x))
                days=(int(start[-2:]),int(pairs[0]['end'][-2:]))
                for dd in range(1,32):
                    subset=values[dd-1,y,x]
                    baseline.append(dict(footprint_start=start,date=f'{year}-05-{dd:02d}',field=field,median=np.nanmedian(subset)))
                    if year==2025 and dd in days:
                        output.append(dict(footprint_start=start,date=f'{year}-05-{dd:02d}',field=field,n_oci=len(subset),
                               n_native_cells=len(set(zip(y,x))),median=np.nanmedian(subset),
                               q25=np.nanquantile(subset,.25),q75=np.nanquantile(subset,.75)))
    for name,rows in [('footprint_physical_context.csv',output),('prior_may_physical_context.csv',baseline)]:
        with (ROOT/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    (ROOT/'footprint_physical_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')

if __name__=='__main__':
    main()
