"""Positional comparisons against published objects and a rolling public TNS broker feed."""
import argparse
import hashlib
import io
import json
import time
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from lxml import html
from campaign import HERE

FIELDS={'cosmos':(150.13,2.35),'ceers':(214.879,52.869),'a370':(39.970,-1.576),
        'whl0137':(24.346,-8.45),'macs0138':(24.515,-21.928),'goods':(53.16,-27.79)}


def download(cache,name,url):
    path=cache/name
    if not path.exists():
        for attempt in range(4):
            try:
                req=Request(url,headers={'User-Agent':'Open-Space-JWST/1.0'})
                with urlopen(req,timeout=60) as response:
                    path.write_bytes(response.read())
                break
            except OSError:
                if attempt==3:raise
                time.sleep(1+attempt)
    data=path.read_bytes()
    return data,dict(source=name,url=url,bytes=len(data),sha256=hashlib.sha256(data).hexdigest())


def sources(cache):
    records=[]; inputs=[]
    data,rec=download(cache,'coulter.html','https://arxiv.org/html/2601.08931')
    inputs.append(rec)
    for row in html.fromstring(data).xpath('//figure[@id="S3.T1"]//tbody/tr'):
        name,ra,dec=[c.text_content().strip() for c in row.xpath('td')[:3]]
        sky=SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
        records.append(['Coulter_COSMOS_Table1',name,sky.ra.deg,sky.dec.deg])
    assert len(records)==68
    url='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'+urlencode({'-source':'J/A+A/714/A22','-out.max':'10000','-out.all':'1'})
    data,rec=download(cache,'dragon.tsv',url);inputs.append(rec)
    tsv=pd.read_csv(io.BytesIO(data),sep='\t',comment='#',dtype=str)
    tsv=tsv[pd.to_numeric(tsv.ID,errors='coerce').notna()]
    for _,row in tsv.iterrows():
        records.append(['Dragon_J_A+A_714_A22','ID '+row.ID.strip(),float(row.RAJ2000),float(row.DEJ2000)])
    assert tsv.ID.nunique()==104
    def get_field(k):
        ra,dec=FIELDS[k]
        url=f'https://transients.fyi/api/sky/markers?ra={ra}&dec={dec}&r=25&window=1Y&source=tns&measure=new&mag=all'
        # Broker has an explicit 1-year rolling window; record coverage for every request.
        data,rec=download(cache,'tns_'+k+'.json',url)
        reply=json.loads(data)
        if reply['truncated'] or reply['source']!='tns':raise ValueError(k)
        return [['TNS_broker_1Y',x['object_id'],x['ra'],x['dec']] for x in reply['objects']],rec
    with ThreadPoolExecutor(5) as pool:
        for rows,rec in pool.map(get_field,FIELDS):
            records.extend(rows);inputs.append(rec)
    kyoto=HERE/'kyoto_field_positions.csv'
    if kyoto.exists():
        records.extend(['Kyoto_TNS_2026',r.name,r.ra,r.dec] for r in pd.read_csv(kyoto).itertuples())
    fink=cache/'fink-fields.csv'
    if not fink.exists() or not (cache/'fink-receipt.json').exists():
        url='https://api.lsst.fink-portal.org/api/v1/resolver'
        query={'resolver':'tns','name_or_id':'','nmax':1000000}
        request=Request(url,data=json.dumps(query).encode(),headers={'Content-Type':'application/json'})
        with urlopen(request,timeout=180) as response:
            data=response.read()
        items=json.loads(data)
        if not (0<len(items)<query['nmax']):raise ValueError('invalid Fink result count')
        centers=[(150.132,2.349),(39.970,-1.576),(24.515,-21.928),(24.346,-8.45),(64.038,-24.068),
                 (101.959,70.248),(268.53,65.2),(53.141,-27.79),(214.879,52.869),(3.595,-30.4)]
        subset=[]
        for obj in items:
            ra,dec=float(obj['f:ra']),float(obj['f:declination'])
            if any(abs(dec-d)<.7 and abs(ra-a)*math.cos(math.radians(d))<.7 for a,d in centers):
                subset.append([obj['f:fullname'],ra,dec])
        pd.DataFrame(subset,columns=['name','ra','dec']).to_csv(fink,index=False)
        receipt=dict(url=url,query=query,sha256=hashlib.sha256(data).hexdigest(),objects=len(items),regional=len(subset),bytes=len(data))
        (cache/'fink-receipt.json').write_text(json.dumps(receipt))
    if fink.exists() and (cache/'fink-receipt.json').exists():
        for row in pd.read_csv(fink).itertuples():
            records.append(['Fink_TNS_resolver',row.name,row.ra,row.dec])
        inputs.append(json.loads((cache/'fink-receipt.json').read_text())|{'source':'Fink_TNS_resolver'})
    return pd.DataFrame(records,columns=['source','identifier','ra','dec']),inputs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache',type=Path);args=parser.parse_args();args.cache.mkdir(parents=True,exist_ok=True)
    refs,inputs=sources(args.cache)
    hits=pd.read_csv(HERE/'coincidences.csv');hitsky=SkyCoord(hits.ra,hits.dec,unit='deg')
    nearest=[]
    for key,group in refs.groupby('source'):
        cat=SkyCoord(group.ra,group.dec,unit='deg');ix,sep,_=hitsky.match_to_catalog_sky(cat)
        for row,j,d in zip(hits.itertuples(),ix,sep.arcsec):
            nearest.append(dict(candidate=row.candidate,catalog=key,nearest=group.iloc[j].identifier,sep_arcsec=d))
    pd.DataFrame(nearest).to_csv(HERE/'nearest_reports.csv',index=False,float_format='%.5f')
    (HERE/'comparison_receipts.json').write_text(json.dumps(inputs,indent=2)+'\n')
    # Provide the field subset for independent re-matching and review.
    refs.drop_duplicates().to_csv(HERE/'comparison_positions.csv',index=False)


if __name__=='__main__':main()
