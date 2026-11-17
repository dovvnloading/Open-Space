"""Measure separations to COSMOS2020 Classic hosts and Gaia DR3 astrometric neighbors."""
import argparse
import io
import json
from urllib.parse import urlencode
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from counterparts import download
from campaign import HERE


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('cache',type=Path)
    cache=parser.parse_args().cache;cache.mkdir(parents=True,exist_ok=True)
    targets=pd.read_csv(HERE/'coincidences.csv')
    config=[
        ('cosmos2020','https://irsa.ipac.caltech.edu/TAP/sync','cosmos2020c','id,ra,dec,lp_type,lp_chis,lp_chi2_best,lp_zpdf,lp_zpdf_l68,lp_zpdf_u68,ez_z500,ez_z160,ez_z840',targets[targets['pair'].str.contains('cosmos')],2),
        ('gaia','https://gea.esac.esa.int/tap-server/tap/sync','gaiadr3.gaia_source','source_id,ra,dec,ref_epoch,pmra,pmdec,parallax,parallax_error',targets,8)]
    summaries=[];receipts=[]
    for name,endpoint,table,cols,hits,radius in config:
        phot=[]
        size = len(hits) if name=='cosmos2020' else 8
        for offset in range(0,len(hits),size):
            chunk=hits.iloc[offset:offset+size]
            constraints=[f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r.ra},{r.dec},{radius/3600:.9f}))" for r in chunk.itertuples()]
            if name=='cosmos2020':
                condition=f"ra BETWEEN {chunk.ra.min()-.001:.6f} AND {chunk.ra.max()+.001:.6f} AND dec BETWEEN {chunk.dec.min()-.001:.6f} AND {chunk.dec.max()+.001:.6f}"
            else:
                condition=' OR '.join(constraints)
            query=f"SELECT TOP 100000 {cols} FROM {table} WHERE {condition}"
            url=endpoint+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query))
            data,receipt=download(cache,f'{name}_b{offset}_v2.csv' if name=='cosmos2020' else f'{name}_b{offset}.csv',url)
            frame=pd.read_csv(io.BytesIO(data),dtype={'source_id':str})
            if 'ra' not in frame:raise ValueError(data[:300])
            phot.append(frame);receipts.append(receipt)
        catalog=pd.concat(phot).drop_duplicates()
        for t in hits.itertuples():
            if len(catalog):
                sky=SkyCoord(catalog.ra,catalog.dec,unit='deg')
                ds=sky.separation(SkyCoord(t.ra,t.dec,unit='deg')).arcsec
                for j in np.argsort(ds):
                    if ds[j]>radius:break
                    summaries.append(dict(candidate=t.candidate,catalog=name,sep_arcsec=ds[j],**catalog.iloc[j].to_dict()))
    pd.DataFrame(summaries).to_csv(HERE/'host_neighbors.csv',index=False,float_format='%.9g')
    (HERE/'host_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')


if __name__=='__main__':main()
