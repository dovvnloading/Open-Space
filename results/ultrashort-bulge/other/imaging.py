"""Check public HST, JWST and Euclid Q2 image footprint metadata at lead positions.

MAST and IRSA are searched in a two-arcsecond aperture. Results are metadata and
can include repeated exposures/products; no PSF contrast limit follows from overlap.
"""
import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from xml.etree import ElementTree as ET
from crossmatch import catalog, CACHE
from fetch import OUT, get, write_csv, SHA


def search(task):
    event,ra,dec,survey=task
    if survey=='MAST':
        req={'service':'Mast.Caom.Filtered.Position','format':'json','params':{
          'columns':'obs_collection,obs_id,instrument_name,filters,t_min,t_max,dataRights,proposal_id',
          'filters':[{'paramName':'obs_collection','values':['HST','JWST']}],
          'position':f'{ra:.8f}, {dec:.8f}, 0.00055556'},'pagesize':500,'page':1}
        url='https://mast.stsci.edu/api/v0/invoke?'+urlencode({'request':json.dumps(req,separators=(',',':'))})
    else:
        url='https://irsa.ipac.caltech.edu/SIA?'+urlencode(dict(COLLECTION='euclid_q2',POS=f'CIRCLE {ra:.8f} {dec:.8f} 0.00055556',MAXREC=500))
    r=dict(event=event,service=survey,url=url)
    path=CACHE/'imaging'/(event+'-'+survey)
    try:
        if not path.exists():path.write_bytes(get(url))
        raw=path.read_bytes();r.update(sha256=SHA(raw),input_bytes=len(raw))
        if survey=='MAST':
            body=json.loads(raw)
            if body['status']!='COMPLETE':
                raise ValueError(body.get('msg',body['status']))
            records=body['data']
            total=body.get('paging',{}).get('rowsFiltered')
            if total is not None and len(records)<total:raise ValueError('incomplete MAST pagination')
            recs=[dict(collection=v['obs_collection'],observation=v['obs_id'],instrument=v['instrument_name'],filter=v['filters'],
                       start_mjd=v['t_min'],stop_mjd=v['t_max'],rights=v['dataRights'],proposal=v['proposal_id']) for v in records]
        else:
            root=ET.fromstring(raw)
            stat=[v for v in root.iter() if v.tag.endswith('INFO') and v.get('name')=='QUERY_STATUS']
            if not stat or any(x.get('value')!='OK' for x in stat):
                raise ValueError(' '.join(x.text or x.get('value','') for x in stat))
            fields=[v.get('name') for v in root.iter() if v.tag.endswith('FIELD')]
            recs=[]
            for line in root.iter():
                if line.tag.endswith('TR'):
                    v=dict(zip(fields,[t.text for t in line]))
                    recs.append(dict(collection='Euclid Q2',observation=v.get('obs_id'),start_mjd=v.get('t_min'),stop_mjd=v.get('t_max'),
                                     product_type=v.get('dataproduct_type'),access_url=v.get('access_url')))
        r.update(status='ok',products=len(recs))
        return r,[dict(event=event,**v) for v in recs]
    except Exception as e:
        r.update(status='error',error=str(e));return r,[]


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('names',nargs='*')
    a=p.parse_args()
    items=catalog();targets=set(a.names)
    if not targets:
        targets={row['event'] for row in csv.DictReader((OUT/'catalog_query_receipts.csv').open())}
    (CACHE/'imaging').mkdir(exist_ok=True)
    tasks=[(n,ra,dec,site) for n,ra,dec,_ in items if n in targets for site in ['MAST','IRSA']]
    with ThreadPoolExecutor(8) as pool: data=list(pool.map(search,tasks))
    write_csv(OUT/'imaging_queries.csv',[r for r,_ in data])
    products=[r for _,ls in data for r in ls]
    if products:write_csv(OUT/'imaging_products.csv',products)
