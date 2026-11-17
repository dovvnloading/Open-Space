"""Compare the repeat-image positions with historical HST stars and COSMOS AGN."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
import argparse
import io
import json
import numpy as np
import pandas as pd
from counterparts import download
from campaign import HERE


def tsv(data,key):
    if key not in data.decode():return pd.DataFrame()
    try:
        a=pd.read_csv(io.BytesIO(data),comment='#',sep='\t',dtype=str)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    if key not in a:return pd.DataFrame()
    return a[pd.to_numeric(a[key],errors='coerce').notna()].copy()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);cache=p.parse_args().cache
    cache.mkdir(parents=True,exist_ok=True);hits=pd.read_csv(HERE/'coincidences.csv')
    hits=hits[hits.pair.str.contains('cosmos|ceers')]
    def hst(h):
        url='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'+urlencode({'-source':'J/ApJS/214/24/3dhstall','-c':f'{h.ra:.7f} {h.dec:+.7f}','-c.rs':40,'-out.max':3000,'-out':'ID,Field,RAJ2000,DEJ2000,S/G,zsp,zpk'})
        data,r=download(cache,'3dhst_v3_'+h.candidate+'.tsv',url)
        r['returned_table']=b'#Table\t' in data
        r['archive_condition']='table_returned' if r['returned_table'] else 'no table returned; lookup cannot exclude a counterpart'
        t=tsv(data,'RAJ2000');rows=[]
        for _,s in t.iterrows():
            ra,dec=float(s.RAJ2000),float(s.DEJ2000)
            sep=np.hypot((h.ra-ra)*np.cos(np.deg2rad(h.dec)),h.dec-dec)*3600
            if sep>2:continue
            rows.append(dict(candidate=h.candidate,catalog='3D-HST',identifier=s.Field.strip()+' '+s.ID.strip(),
                             sep_arcsec=sep,ra=ra,dec=dec,star_flag=s['S/G'],z_spec=s.zsp.strip(),z_peak=s.zpk.strip()))
        return rows,r
    with ThreadPoolExecutor(6) as pool:batches=list(pool.map(hst,hits.itertuples()))
    rows=[row for data,_ in batches for row in data];receipts=[r for _,r in batches]
    u='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'+urlencode({'-source':'J/A+A/664/A117','-out.max':'10000','-out.add':'_RAJ2000,_DEJ2000','-out.all':'1'})
    dat,receipt=download(cache,'vst_all.tsv',u);receipts.append(receipt)
    vst=tsv(dat,'RAJ2000').rename(columns={'RAJ2000':'ra','DEJ2000':'dec','ID':'name'});assert len(vst)==677
    q='SELECT TOP 1000 * FROM cosmos_chandra_2016 WHERE ra BETWEEN 150.07 AND 150.19 AND dec BETWEEN 2.23 AND 2.43'
    u='https://irsa.ipac.caltech.edu/TAP/sync?'+urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':q})
    data,receipt=download(cache,'xray_irsa.csv',u);receipts.append(receipt)
    xray=pd.read_csv(io.BytesIO(data));xray=xray.rename(columns={xray.columns[0]:'name'})
    for cat,t in [('VST-COSMOS_2022',vst),('Chandra_COSMOS_Legacy',xray)]:
        for h in hits[hits.pair.str.contains('cosmos')].itertuples():
            sep=3600*np.hypot((t.ra.astype(float)-h.ra)*np.cos(np.deg2rad(h.dec)),t.dec.astype(float)-h.dec)
            j=sep.idxmin()
            rows.append(dict(candidate=h.candidate,catalog=cat,identifier=t.loc[j,'name'],sep_arcsec=sep[j],ra=t.ra[j],dec=t.dec[j]))
    pd.DataFrame(rows).to_csv(HERE/'archival_neighbors.csv',index=False,float_format='%.9g')
    (HERE/'archival_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')


if __name__=='__main__':main()
