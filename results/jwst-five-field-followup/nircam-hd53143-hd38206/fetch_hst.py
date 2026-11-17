import argparse,csv,json,hashlib,shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
from datetime import datetime,timezone
def main(out,cache):
    d=json.loads((out/'HD38206-hst.json').read_text())
    rows=[r for r in d['records'] if r['obs_id'] in ('j9f919030','j9f920020','j9f921030')]
    assert len(rows)==3
    def get(r):
        u='https://mast.stsci.edu/api/v0.1/Download/file?uri='+r['dataURL'];p=cache/r['dataURL'].split('/')[-1]
        with urlopen(u,timeout=150) as h,p.with_suffix('.part').open('wb') as o:shutil.copyfileobj(h,o)
        p.with_suffix('.part').replace(p)
        with p.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
        return {**r,'url':u,'bytes':p.stat().st_size,'sha256':sha,'utc':datetime.now(timezone.utc).isoformat()}
    with ThreadPoolExecutor(3) as pool:rows=list(pool.map(get,rows))
    (out/'hst_inputs.json').write_text(json.dumps(rows,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);p.add_argument('cache',type=Path);a=p.parse_args();main(a.out,a.cache)
