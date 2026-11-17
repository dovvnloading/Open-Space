"""Get per-roll clear-aperture mosaics for the June 2024 independent epoch of HD 53143."""
import csv, io, json, hashlib, argparse, shutil
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor

def main(cache,out):
    q=dict(service='Mast.Caom.Filtered',format='json',params={'columns':'obsid,obs_id,target_name,filters,calib_level,dataRights,t_min,t_max,t_obs_release,dataURL','filters':[{'paramName':'obs_collection','values':['JWST']},{'paramName':'proposal_id','values':['3973']},{'paramName':'calib_level','values':[3]}]},pagesize=1000,page=1)
    url='https://mast.stsci.edu/api/v0/invoke?'+urlencode({'request':json.dumps(q)})
    with urlopen(url,timeout=90) as r:body=r.read()
    data=json.loads(body);assert data['status']=='COMPLETE'
    (out/'wide_CAOM.json').write_text(json.dumps(dict(request=q,response_sha256=hashlib.sha256(body).hexdigest(),retrieved_utc=datetime.now(timezone.utc).isoformat(),rows=data['data']),indent=2)+'\n')
    rows=[r for r in data['data'] if 'F444W' in str(r['filters']) and r['obs_id'].startswith(('jw03973-o006','jw03973-o007')) and r['dataURL'] and r['dataURL'].endswith('_i2d.fits')]
    def fetch(r):
        u='https://mast.stsci.edu/api/v0.1/Download/file?uri='+r['dataURL'];p=cache/r['dataURL'].split('/')[-1]
        if not p.exists():
            with urlopen(u,timeout=240) as h,p.with_suffix('.part').open('wb') as fo:
                shutil.copyfileobj(h,fo)
            p.with_suffix('.part').replace(p)
        with p.open('rb') as f: sha=hashlib.file_digest(f,'sha256').hexdigest()
        return {**r,'url':u,'bytes':p.stat().st_size,'sha256':sha}
    with ThreadPoolExecutor(2) as pool:inputs=list(pool.map(fetch,rows))
    (out/'wide_inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('cache',type=Path);p.add_argument('out',type=Path);a=p.parse_args();main(a.cache,a.out)
