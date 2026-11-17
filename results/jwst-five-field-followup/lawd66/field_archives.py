"""Refresh/replay 45-arcsecond cones on the LAWD 66 field through the original archive helpers."""
import argparse
import concurrent.futures
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import urlopen

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('repo',type=Path)
p.add_argument('--refresh',action='store_true')
args=p.parse_args()
out=Path(__file__).resolve().parent/'archives'
out.mkdir(exist_ok=True)
sys.path.insert(0,str(args.repo/'results/jwst-nearby-companions/miri'))
from catalogs import get_service, SERVICES
ra,dec=257.86320,-14.79990
meta={'epochs':[{'ra':ra,'dec':dec}]}
q="SELECT obs_id,obs_collection,instrument_name,dataproduct_type,t_min,t_max FROM ivoa.obscore WHERE 1=INTERSECTS(s_region,CIRCLE('ICRS',257.86320,-14.79990,0.0125))"
urls={'HST':('csv','https://hst.esac.esa.int/tap-server/tap/sync?'+urlencode(dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",QUERY=q))),'PS1':('csv','https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv?'+urlencode(dict(ra=ra,dec=dec,radius=0.0125,columns='[objID,raMean,decMean,epochMean,gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,yMeanPSFMag,nDetections]')))}
def fetch(n):
    if n in SERVICES:
        get_service(('LAWD66',n,meta,out,args.refresh))
        return
    ext,url=urls[n]
    dst=out/f'{n}_response.{ext}'
    stamp=out/f'{n}_receipt.json'
    if args.refresh:
        with urlopen(url,timeout=100) as response:
            raw=response.read()
        if n=='HST' and not raw.startswith(b'obs_id,'): raise ValueError(raw[:120])
        if n=='PS1' and not raw.startswith(b'objID,'): raise ValueError(raw[:120])
        dst.write_bytes(raw)
        stamp.write_text(json.dumps(dict(url=url,utc=datetime.now(timezone.utc).isoformat(),sha256=hashlib.sha256(raw).hexdigest()),indent=2)+'\n')
    if hashlib.sha256(dst.read_bytes()).hexdigest()!=json.loads(stamp.read_text())['sha256']: raise ValueError(dst)
with concurrent.futures.ThreadPoolExecutor(5) as pool:
    futures={pool.submit(fetch,n):n for n in list(SERVICES)+list(urls)}
    for job in concurrent.futures.as_completed(futures):
        try: job.result();print(futures[job],'success',flush=True)
        except Exception as e: print(futures[job],repr(e),flush=True); raise
