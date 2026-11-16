"""Fetch daily NOAA physical context at the nearest grid cells. Python 3, standard library."""
import csv
from concurrent.futures import ThreadPoolExecutor
from datetime import date, timedelta
import hashlib
import json
import math
from pathlib import Path
import re
from urllib.parse import quote
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
PSL = 'https://psl.noaa.gov/thredds/dodsC/Datasets/'
GEO = 'https://cwcgom.aoml.noaa.gov/thredds/dodsC/OCEAN_GEOSTROPHIC_CURRENTS/CURRENTS.nc'
SITES = {
    'usc_seaprism': (33.56371, -118.11782, [('2024-03-11','2024-03-16'),('2024-04-10','2024-04-15')]),
    'florida_atlantic': (29.0, -79.84, [('2025-05-13','2025-05-21')]),
}
GAUSSIAN_LATITUDES = {29: 33.3328, 30:31.4281,31:29.5234,32:27.6186}
FIELDS = ['anom','u_current','v_current','uwnd','vwnd']

def prepare(site, first, last, field):
    latitude, longitude, _ = SITES[site]
    if field == 'anom':
        index0 = date(first.year, 1, 1)
        origin, scale = date(1800, 1, 1), 1.
        base = PSL+f'noaa.oisst.v2.highres/sst.day.anom.{first.year}.nc'
        # A geographic tie uses the lower latitude.
        y = math.ceil((latitude+89.875)/.25-.5)
        x = round(((longitude%360)-.125)/.25)
    elif field.endswith('current'):
        index0 = date(2016, 7, 1)
        origin, scale = date(1970, 1, 1), 86400.
        base = GEO
        y = math.ceil((latitude+64.7)/.2-.5)
        x = round((longitude+179.9)/.2)
    else:
        index0 = date(first.year, 1, 1)
        origin, scale = date(1800, 1, 1), 24.
        base = PSL+f'ncep.reanalysis.dailyavgs/surface_gauss/{field}.10m.gauss.{first.year}.nc'
        y = min(GAUSSIAN_LATITUDES,key=lambda i:abs(GAUSSIAN_LATITUDES[i]-latitude))
        x = round((longitude%360)/1.875)
    constraint = f'{field}[{(first-index0).days}:1:{(last-index0).days}][{y}][{x}]'
    return base+'.ascii?'+quote(constraint), origin, scale

def fetch(task):
    site,first,last,field=task
    url,origin,scale=prepare(*task)
    with urlopen(url,timeout=75) as res: raw=res.read()
    payload=raw.decode().split('---------------------------------------------\n')[1]
    def coord(key):
        match=re.search(rf'{field}\.{key}\[\d+\]\n([^\n]+)',payload)
        return [float(t) for t in match.group(1).split(',')]
    data = re.findall(r'^\[\d+\]\[0\], ([^\n]+)',payload,re.M)
    lat,lon=coord('lat')[0],(coord('lon')[0]+180)%360-180
    rows=[]
    for i,(t,value) in enumerate(zip(coord('time'),data,strict=True)):
        target=first+timedelta(days=i)
        if t != (target-origin).days*scale: raise ValueError('Timestamp mismatch')
        rows.append(dict(site=site,date=target.isoformat(),field=field,grid_lat=round(lat,6),
                         grid_lon=round(lon,6),value=round(float(value),6)))
    if len(rows)!=(last-first).days+1: raise ValueError('Incomplete range')
    return rows,dict(site=site,first=str(first),last=str(last),field=field,
                     url=url,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))

def main():
    tasks=[(site,date.fromisoformat(a),date.fromisoformat(b),field)
           for site,(_,_,windows) in SITES.items() for a,b in windows for field in FIELDS]
    with ThreadPoolExecutor(4) as pool: results=list(pool.map(fetch,tasks))
    rows=[row for subset,_ in results for row in subset]
    with (ROOT/'daily_physical_context.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    (ROOT/'physical_receipts.json').write_text(json.dumps([r for _,r in results],indent=2)+'\n')

if __name__=='__main__':
    main()
