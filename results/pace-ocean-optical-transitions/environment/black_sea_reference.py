"""Count matched-duration March changes in the Section-7 multi-year reanalysis."""
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import date, timedelta
import json

import numpy as np
from regional_hydrography import BLOSC, read, RECEIPTS, ROOT

PRODUCT='BLKSEA_MULTIYEAR_PHY_007_004'
POINT=(44.5458,29.4466)

def series(field):
    family={'so':'sal','thetao':'temp'}[field]
    ds=f'cmems_mod_blk_phy-{family}_my_2.5km_P1D-m_202411'
    stac=f'https://stac.marine.copernicus.eu/metadata/{PRODUCT}/{ds}/dataset.stac.json'
    base=json.loads(read(stac))['assets']['geoChunked']['href']
    meta=json.loads(read(base+'/.zmetadata'))['metadata']
    def load(v,indices):
        spec=meta[v+'/.zarray'];raw=read(base+'/'+v+'/'+'.'.join(map(str,indices)))
        out=np.empty(spec['chunks'],dtype=spec['dtype'])
        assert spec['compressor']['id']=='blosc'
        assert BLOSC.blosc_decompress(raw,out.ctypes.data,out.nbytes)==out.nbytes
        return out
    def axis(v):
        sp=meta[v+'/.zarray'];n=sp['shape'][0];c=sp['chunks'][0]
        return np.concatenate([load(v,(i,)) for i in range((n+c-1)//c)])[:n]
    sy,sx=POINT
    lat,lon,elev,t=[axis(v) for v in ('latitude','longitude','elevation','time')]
    y=int(np.abs(lat-sy).argmin());x=int(np.abs(lon-sx).argmin());z=int(np.abs(elev).argmin())
    spec=meta[field+'/.zarray']
    ct,cz,cy,cx=spec['chunks']
    def timedata(i):
        return load(field,(i,z//cz,y//cy,x//cx))[:,z%cz,y%cy,x%cx]
    with ThreadPoolExecutor(3) as pool:
        values=np.concatenate(list(pool.map(timedata,range((len(t)+ct-1)//ct))))[:len(t)]
    assert meta['time/.zattrs']['units']=='seconds since 1970-01-01'
    dates=[date(1970,1,1)+timedelta(days=int(ti)//86400) for ti in t]
    assert all((d1-d0).days==1 for d0,d1 in zip(dates,dates[1:]))
    return {d:float(v) for d,v in zip(dates,values) if d.month==3 and d.year<=2025 and np.isfinite(v) and v<1e10}

def main():
    with ThreadPoolExecutor(2) as pool: sst,sal=pool.map(series,['thetao','so'])
    dts=sorted(set(sst)&set(sal))
    output=[]
    for d in dts:
        after=d+timedelta(days=3)
        if after not in sal or after not in sst: continue
        output.append(dict(start=d.isoformat(),end=after.isoformat(),
                           change_temperature_c=round(sst[after]-sst[d],6),
                           change_salinity=round(sal[after]-sal[d],6)))
    with (ROOT/'section7_march_physical_reference.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=output[0]);w.writeheader();w.writerows(output)
    (ROOT/'section7_march_reference_receipts.json').write_text(json.dumps([RECEIPTS[k] for k in sorted(RECEIPTS)],indent=2)+'\n')

if __name__=='__main__': main()
