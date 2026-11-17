"""Retrieve and hash the public SHA IRAC2 mosaics and associated coverage maps."""
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE=Path(__file__).resolve().parent
SIA='https://irsa.ipac.caltech.edu/SIA'
PARAMS={'POS':'CIRCLE 335.7815 -17.6125 0.01','BAND':'0.0000040 0.0000051',
        'COLLECTION':'spitzer_sha','FACILITY':'Spitzer Space Telescope','RESPONSEFORMAT':'CSV'}
OBSLOG='https://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/spitzer_obslog.txt'

def digest(blob):return hashlib.sha256(blob).hexdigest()

def main(cache):
    cache.mkdir(exist_ok=True,parents=True)
    with urlopen(SIA+'?'+urlencode(PARAMS),timeout=120) as f:body=f.read()
    receipt={'utc':datetime.now(timezone.utc).isoformat(),'endpoint':SIA,'params':PARAMS,'sha256_response':digest(body)}
    rows=list(csv.DictReader(io.StringIO(body.decode())))
    receipt['rows']=len(rows)
    select={r['access_url'] for r in rows if r['instrument_name']=='IRAC' and r['energy_bandpassname']=='IRAC2' and
            '/pbcd/' in r['access_url'] and r['access_url'].endswith(('_maic.fits','_mcov.fits')) and
            any('/r'+aid+'/' in r['access_url'] for aid in ('57886720','65885696'))}
    if len(select)!=4:raise ValueError(f'Expected four PBCD inputs; archive returned {len(select)}')
    inputs=[]
    for url in sorted(select):
        dest=cache/url.rsplit('/',1)[-1]
        if not dest.exists():
            with urlopen(url,timeout=120) as f:body=f.read()
            if not body.startswith(b'SIMPLE  ='): raise ValueError('FITS signature missing')
            dest.write_bytes(body)
        body=dest.read_bytes()
        inputs.append({'url':url,'filename':dest.name,'bytes':len(body),'sha256':digest(body)})
    with urlopen(OBSLOG,timeout=90) as f:log=f.read()
    receipt['executed_log']={'url':OBSLOG,'sha256':digest(log),'matching_lines':[line.strip() for line in log.decode().splitlines() if '57886720' in line or '65885696' in line]}
    receipt['files']=inputs
    (HERE/'spitzer-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path);main(ap.parse_args().cache)
