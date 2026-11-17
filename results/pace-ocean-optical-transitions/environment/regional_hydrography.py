"""Read versioned Copernicus 4.2/2.5-km analyses at the AAOT and Section-7 stations."""
import csv
import ctypes
import ctypes.util
from datetime import date, timedelta
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from time import sleep
from urllib.request import urlopen

import numpy as np

from physical_context import SITES, fetch

ROOT=Path(__file__).resolve().parent
OCEANS={
    'AAOT':('med','MEDSEA_ANALYSISFORECAST_PHY_006_013','4.2km',45.3139,12.5083,date(2025,11,4),date(2025,11,9)),
    'Section-7':('blk','BLKSEA_ANALYSISFORECAST_PHY_007_001','2.5km',44.5458,29.4466,date(2026,3,7),date(2026,3,15)),
}
RECEIPTS={}
BLOSC=ctypes.CDLL(ctypes.util.find_library('blosc') or 'libblosc.so.1')
BLOSC.blosc_decompress.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_size_t]

def read(url):
    cache=ROOT/'cache'
    cache.mkdir(exist_ok=True)
    path=cache/hashlib.sha256(url.encode()).hexdigest()
    if path.exists(): data=path.read_bytes()
    else:
        for attempt in range(4):
            try:
                with urlopen(url,timeout=50) as r: data=r.read()
                break
            except OSError:
                if attempt==3: raise
                sleep(2**attempt)
        path.write_bytes(data)
    RECEIPTS[url]=dict(url=url,bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
    return data

def fields(site,family,variables):
    tag,product,res,slat,slon,start,end=OCEANS[site]
    if tag=='blk' and family=='tem': family='temp'
    dataset=f'cmems_mod_{tag}_phy-{family}_anfc_{res}_P1D-m_202511'
    stac=f'https://stac.marine.copernicus.eu/metadata/{product}/{dataset}/dataset.stac.json'
    base=json.loads(read(stac))['assets']['geoChunked']['href']
    meta=json.loads(read(base+'/.zmetadata'))['metadata']
    @lru_cache(maxsize=6)
    def chunk(variable,indices):
        spec=meta[variable+'/.zarray']
        path=variable+'/'+spec.get('dimension_separator','.').join(map(str,indices))
        raw=read(base+'/'+path)
        out=np.empty(spec['chunks'],dtype=spec['dtype'])
        assert spec['compressor']['id']=='blosc' and spec['order']=='C'
        size=BLOSC.blosc_decompress(raw,out.ctypes.data,out.nbytes)
        if size!=out.nbytes: raise ValueError('Decompression failure')
        return out
    def coordinate(v):
        spec=meta[v+'/.zarray']
        n=spec['shape'][0];ch=spec['chunks'][0]
        return np.concatenate([chunk(v,(i,)) for i in range((n+ch-1)//ch)])[:n]
    lat,lon,depth,times=[coordinate(v) for v in ('latitude','longitude','elevation','time')]
    y=int(np.abs(lat-slat).argmin());x=int(np.abs(lon-slon).argmin());z=int(np.abs(depth).argmin())
    unit,epoch=meta['time/.zattrs']['units'].split(' since ')
    scale={'seconds':86400,'minutes':1440}[unit]
    dates=[date.fromisoformat(epoch)+timedelta(days=int(t)//scale) for t in times]
    positions=np.flatnonzero([(start<=d<=end) for d in dates])
    if len(positions)!=(end-start).days+1: raise ValueError('Incomplete coverage')
    rows=[]
    for i in positions:
        for var in variables:
            spec=meta[var+'/.zarray']
            def value(j,k):
                index=np.array([i,z,j,k]);c=np.array(spec['chunks'])
                v=float(chunk(var,tuple(index//c))[tuple(index%c)])
                return np.nan if v==spec['fill_value'] else v
            center=value(y,x)
            neighbors=[value(j,k) for j in range(y-2,y+3) for k in range(x-2,x+3)]
            rows.append(dict(site=site,date=dates[i].isoformat(),field=var,depth_m=round(-float(depth[z]),4),
                lat=round(float(lat[y]),6),lon=round(float(lon[x]),6),cell=round(center,6),
                n_water_5x5=int(np.isfinite(neighbors).sum()),median_5x5=round(float(np.nanmedian(neighbors)),6),
                q25_5x5=round(float(np.nanquantile(neighbors,.25)),6),q75_5x5=round(float(np.nanquantile(neighbors,.75)),6)))
    return rows

def main():
    from concurrent.futures import ThreadPoolExecutor
    jobs=[(site,fam,vars) for site in OCEANS for fam,vars in [('cur',['uo','vo']),('sal',['so']),('tem',['thetao'])]]
    with ThreadPoolExecutor(6) as pool:
        rows=[r for sub in pool.map(lambda j:fields(*j),jobs) for r in sub]
    for site,(_,_,_,lat,lon,start,end) in OCEANS.items():
        SITES[site]=(lat,lon,[])
        point,receipt=fetch((site,start,end,'anom'))
        RECEIPTS[receipt['url']]=receipt
        for p in point: rows.append(dict(site=site,date=p['date'],field='oisst_anomaly',depth_m=0,
            lat=p['grid_lat'],lon=p['grid_lon'],cell=p['value'],n_water_5x5='',median_5x5='',q25_5x5='',q75_5x5=''))
    with (ROOT/'regional_hydrography.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    (ROOT/'regional_physics_receipts.json').write_text(json.dumps([RECEIPTS[k] for k in sorted(RECEIPTS)],indent=2)+'\n')

if __name__=='__main__': main()
