from urllib.request import urlopen
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import csv,sys
D=Path(__file__).resolve().parent
obs={'a2744':(['jw03073-o008_t015','jw06434-c1021_t000'],['f200w','f444w']),'m0647':(['jw01433-o010_t001','jw01433-o020_t001'],['f200w'])}
def fetch(arg):
    field,obs,band=arg; root=obs+'_nircam_clear-'+band
    p,o=obs.split('_')[0].split('-')
    url=f'https://stpubdata.s3.amazonaws.com/jwst/public/{p}/L3/t/{o}/{root}_cat.ecsv'
    b=urlopen(url,timeout=30).read()
    (D/'catalog-cache'/(root+'.ecsv')).write_bytes(b)
    return dict(field=field,product=root,url=url,bytes=len(b),sha256=sha256(b).hexdigest())
if __name__=='__main__':
 (D/'catalog-cache').mkdir(exist_ok=True)
 z=list(ThreadPoolExecutor(4).map(fetch,[(f,o,b) for f,(oo,bb) in obs.items() for o in oo for b in bb]))
 with (D/'catalog_inputs.csv').open('w') as h:
  w=csv.DictWriter(h,fieldnames=z[0]);w.writeheader();w.writerows(z)
