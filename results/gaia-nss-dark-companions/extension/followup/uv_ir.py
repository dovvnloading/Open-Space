"""Check GALEX and WISE measurements at the giant's Gaia position."""
import csv
import hashlib
import json
import urllib.parse
import urllib.request
from pathlib import Path

HERE=Path(__file__).resolve().parent
s=next(csv.DictReader((HERE/'source.csv').open()))
ra,dec=map(float,(s['ra'],s['dec']))
request=dict(service='Mast.Galex.Catalog',params=dict(ra=ra,dec=dec,radius=3/3600),format='json')
where=f" WHERE 1=CONTAINS(POINT('J2000',ra,dec),CIRCLE('J2000',{ra},{dec},0.00083333333))"
q="SELECT designation,ra,dec,w1mpro,w1sigmpro,w2mpro,w2sigmpro,w3mpro,w3sigmpro,w4mpro,ph_qual,cc_flags FROM allwise_p3as_psd"+where
urls=dict(galex='https://mast.stsci.edu/api/v0/invoke?'+urllib.parse.urlencode({'request':json.dumps(request)}),
          allwise='https://irsa.ipac.caltech.edu/TAP/sync?'+urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q)))
receipts=[]
for key,url in urls.items():
    with urllib.request.urlopen(url,timeout=70) as r: data=r.read()
    if key=='galex':
        response=json.loads(data)
        if response.get('status')!='COMPLETE':raise ValueError('GALEX response did not complete')
        rows=len(response['data'])
    else:
        rows=len(list(csv.DictReader(data.decode().splitlines())))
    filename=key+('.json' if key=='galex' else '.csv')
    (HERE/filename).write_bytes(data)
    receipts.append(dict(file=filename,url=url,rows=rows,sha256=hashlib.sha256(data).hexdigest(),radius_arcsec=3))
(HERE/'uv_ir_inputs.json').write_text(json.dumps(receipts,indent=2)+'\n')
