"""Public Rubin alert retrieval following Open-Space's ALeRCE TAP campaign."""
import argparse, hashlib, json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone

TAP='https://tap.alerce.online/tap/sync?'
FINK='https://api.lsst.fink-portal.org/api/v1/conesearch?'
queries={
'schema':"SELECT column_name FROM TAP_SCHEMA.columns WHERE table_name='alerce_tap.detection'",
'objects':"SELECT TOP 30000 o.oid,o.n_det,o.meanra,o.meandec,o.firstmjd,o.lastmjd FROM alerce_tap.object AS o WHERE o.sid=1 AND o.meanra BETWEEN 34.64 AND 35.43 AND o.meandec BETWEEN -49.96 AND -49.22",
'alerts':"SELECT TOP 100000 d.oid,d.measurement_id,d.mjd,d.ra,d.dec,d.band FROM alerce_tap.detection AS d JOIN alerce_tap.object AS o ON d.sid=o.sid AND d.oid=o.oid WHERE o.sid=1 AND o.meanra BETWEEN 34.64 AND 35.43 AND o.meandec BETWEEN -49.96 AND -49.22 AND d.mjd BETWEEN 60900 AND 61362"}

def run(name,url,dest):
    receipt=dict(name=name,url=url,utc=datetime.now(timezone.utc).isoformat())
    try:
        with urlopen(url,timeout=150) as f:b=f.read();receipt['status']=f.status
        (dest/(name+'.data')).write_bytes(b);receipt['sha256']=hashlib.sha256(b).hexdigest();receipt['bytes']=len(b);print(name,len(b),repr(b[:230]),flush=True)
    except Exception as e:receipt['error']=repr(e);print(name,repr(e),flush=True)
    (dest/(name+'.json')).write_text(json.dumps(receipt,indent=1)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('dest',type=Path);a=p.parse_args();a.dest.mkdir(exist_ok=True,parents=True)
    urls={n:TAP+urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','MAXREC':100000,'QUERY':q}) for n,q in queries.items()}
    urls['fink']=FINK+urlencode({'ra':35.045,'dec':-49.62,'radius':1700,'columns':'r:diaObjectId,r:midpointMjdTai,r:ra,r:dec,r:psfFlux,r:psfFluxErr','output-format':'csv'})
    with ThreadPoolExecutor(4) as ex:list(ex.map(lambda kv:run(*kv,a.dest),urls.items()))
