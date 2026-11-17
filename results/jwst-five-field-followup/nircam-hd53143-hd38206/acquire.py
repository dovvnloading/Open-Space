"""Materialize the public NIRCam exposures identified in the merged inventory.
Run with --parent PATH to the merged jwst-nearby-companions results.
Known calibration inputs are authenticated by their original manifest.
"""
import argparse, csv, gzip, hashlib, io, json, shutil, sys, time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

def main(parent,out,cache):
    sys.path.insert(0,str(parent/'nircam'))
    from nircam_fetch import digest
    cache.mkdir(exist_ok=True,parents=True)
    out.mkdir(exist_ok=True,parents=True)
    with gzip.open(parent/'inventory/archive/coron.csv.gz','rt') as f:
        rows=[r for r in csv.DictReader(f) if '-' not in r['fileSetName'] and (r['program'],r['observtn']) in
              {('3973',i) for i in ['1','3','8']}|{('6012',i) for i in ['141','142','150','151']}]
    # The fileset inventory has both combined optical elements on a single row.
    base='https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/'
    items=[]
    for r in sorted(rows,key=lambda x:x['fileSetName']):
        assert r['access'].lower() in ('public','1'),r['access']
        for detector in ('nrcalong','nrca2'):
            name=r['fileSetName']+'_'+detector+'_calints.fits'
            items.append(dict(filename=name,url=base+name,role='reference' if r['is_psf']=='t' else 'science',
                         program=r['program'],obs=r['observtn'],target=r['targname'],release=r['publicReleaseDate'],
                         fileset_id=r['ArchiveFileID'],expected_sha256=''))
    with (parent/'nircam/mcc124_inputs.csv').open() as f:
        for r in csv.DictReader(f):
            if r['filename'].startswith('jw05835'):continue
            items.append(dict(filename=r['filename'],url=r['url'],role=r['role'],program='',obs='',target=r['target'],
                              release=r['release_mjd'],fileset_id='',expected_sha256=r['sha256']))
    def download(rec):
        path=cache/rec['filename']
        params={}
        if not path.exists():
            for attempt in range(4):
                try:
                    with urlopen(rec['url'],timeout=180) as res,path.with_suffix('.part').open('wb') as fo:
                        shutil.copyfileobj(res,fo,1024*1024)
                        params=dict(http_status=res.status,last_modified=res.headers.get('last-modified',''),content_length=res.headers.get('content-length',''),retrieved_utc=datetime.now(timezone.utc).isoformat())
                    partial=path.with_suffix('.part')
                    if params['content_length']:assert partial.stat().st_size==int(params['content_length'])
                    partial.replace(path)
                    (cache/(path.name+'.receipt')).write_text(json.dumps(params))
                    break
                except OSError:
                    if attempt==3:raise
                    time.sleep(2**attempt)
        params=json.loads((cache/(path.name+'.receipt')).read_text())
        sha=digest(path)
        if rec['expected_sha256'] and sha!=rec['expected_sha256']: raise ValueError(f'{path} differs from original hash: {sha}')
        return dict(**rec,bytes=path.stat().st_size,sha256=sha,**params)
    with ThreadPoolExecutor(8) as pool:
        futures=[pool.submit(download,i) for i in items]
        ok=[]; errors=[]
        for i,f in zip(items,futures):
            try: ok.append(f.result())
            except Exception as err:errors.append(dict(filename=i['filename'],error=str(err)))
            if ok:
                with (out/'inputs.csv').open('w') as f:
                    w=csv.DictWriter(f,fieldnames=list(ok[0]));w.writeheader();w.writerows(ok)
            (out/'download-errors.json').write_text(json.dumps(errors,indent=2)+'\n')
    if errors:raise RuntimeError(errors)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cache',type=Path,required=True)
    args=p.parse_args();main(args.parent,args.out,args.cache)
