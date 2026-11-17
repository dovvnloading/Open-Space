"""Fetch actual W1 science/masks in the same final-mission frames selected for W2."""
import sys,io,csv,hashlib,gzip,json
from pathlib import Path
from urllib.parse import urlencode
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from astropy.io import fits
from acquire import HERE,load_url

def main(cache):
    w2=json.loads((HERE/'receipts'/'Phan.json').read_text());folder=cache/'Phan-W1';folder.mkdir(parents=True,exist_ok=True)
    u='https://irsa.ipac.caltech.edu/SIA?'+urlencode(dict(COLLECTION='neowiser',TIME='59945 60532',POS='CIRCLE 35.05000 -49.63000 0.04',BAND='3.4e-6',RESPONSEFORMAT='CSV',MAXREC=30000));raw=load_url(u);(folder/'query.csv').write_bytes(raw)
    rows=list(csv.DictReader(io.StringIO(raw.decode())));science={x['obs_id']:x for x in rows if x.get('energy_bandpassname')=='W1' and '-w1-int-' in x['access_url']};assert science
    led=dict(target='Phan-W1',band='W1',query=u,query_sha256=hashlib.sha256(raw).hexdigest(),rows=len(rows),science_count=len(science),center=w2['center'],retrieved_utc=datetime.now(timezone.utc).isoformat(),products=[])
    for q in w2['products']:
        if q['id'] in science:
            x=science[q['id']];led['products'].append(dict(id=q['id'],visit=q['visit'],mjd=(float(x['t_min'])+float(x['t_max']))/2,url=x['access_url']))
    def fetch(p):
        for k in ('int','msk'):
            suffix=p['url'].split('/p1bm_frm/')[1].replace('-int-','-'+k+'-')+('.gz' if k=='msk' else '');url='https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/neowiser/images/p1bm_frm/'+suffix
            p[k+'_url']=url;p[k+'_original_url']=url
            try:
                b=load_url(url);d=gzip.decompress(b) if b[:2]==bytes([31,139]) else b;assert d[:6]==b'SIMPLE';path=folder/(p['id']+'_'+k+'.fits');path.write_bytes(d);h=fits.getheader(path);p[k+'_sha256']=hashlib.sha256(d).hexdigest();p[k+'_original_sha256']=hashlib.sha256(b).hexdigest();p[k+'_shape']=[h['NAXIS2'],h['NAXIS1']]
            except Exception as e:p[k+'_error']=repr(e)
    with ThreadPoolExecutor(5) as ex:list(ex.map(fetch,led['products']))
    (HERE/'receipts'/'Phan-W1.json').write_text(json.dumps(led,indent=1)+'\n')
if __name__=='__main__':main(Path(sys.argv[1]))
