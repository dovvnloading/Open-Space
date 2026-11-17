"""Resolve the observation intervals of the products used in the catalog screen."""
import argparse
import csv
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from astropy.time import Time
import pandas as pd
from campaign import HERE
from counterparts import download

CENTERS={'cosmos':(150.132,2.349),'ceers':(214.8788,52.8688),'a370':(39.970,-1.576),
         'macs0138':(24.515,-21.928),'whl0137':(24.346,-8.45)}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path)
    cache=p.parse_args().cache;cache.mkdir(parents=True,exist_ok=True)
    selected=set(pd.read_csv(HERE/'catalog_inputs.csv')['product'])
    def get(item):
        field,(ra,dec)=item
        query=f"SELECT TOP 10000 obs_id,t_min,t_max,access_url FROM ivoa.ObsCore WHERE calib_level=3 AND instrument_name='NIRCAM/IMAGE' AND dataproduct_type='image' AND 1=CONTAINS(POINT('ICRS',{ra},{dec}),s_region) ORDER BY t_max DESC"
        url='https://jwst.esac.esa.int/server/tap/sync?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query))
        data,receipt=download(cache,'inventory_'+field+'.csv',url)
        rows=[]
        for row in csv.DictReader(io.StringIO(data.decode())):
            if row['obs_id'] in selected:
                t0,t1=float(row['t_min']),float(row['t_max'])
                rows.append(dict(field=field,product=row['obs_id'],mjd_start=t0,mjd_end=t1,utc_start=Time(t0,format='mjd').isot,utc_end=Time(t1,format='mjd').isot))
        return rows,receipt
    with ThreadPoolExecutor(5) as pool:batches=list(pool.map(get,CENTERS.items()))
    out=pd.DataFrame([r for rows,_ in batches for r in rows]).drop_duplicates(subset='product').sort_values('product')
    if set(out['product'])!=selected:raise ValueError('missing product dates')
    out.to_csv(HERE/'observing_intervals.csv',index=False)
    (HERE/'inventory_receipts.json').write_text(json.dumps([r for _,r in batches],indent=2)+'\n')


if __name__=='__main__':main()
