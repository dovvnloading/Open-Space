"""Fetch the northern ecliptic pole selection, counterpart inputs and DR1 tile.

Run in the repository environment with Python 3, numpy, astropy. Retain the
small inputs and receipts. The FITS cube is downloaded to localcache/.
"""
import csv
import hashlib
import json
import os
import sys
import time
import tempfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
FIELD = os.environ.get('SPHEREX_FIELD','nep')
if FIELD not in ('nep','wide'): raise ValueError(FIELD)
ROOT = HERE if FIELD=='nep' else HERE/'wide'
REPO = next(p for p in [HERE.parent.parent, HERE.parent/'repo'] if (p/'results/neowise-nuclear-flares/spherex_tools.py').exists())
sys.path.insert(0, str(REPO/'results/neowise-nuclear-flares'))
from spherex_tools import tile_url, healpix_nest

RA, DEC = (270.,66.56) if FIELD=='nep' else (209.53125,30.)
TILE = healpix_nest(5, RA, DEC)
CENTER = f'{RA}, {DEC}'
WISE = f"""SELECT designation,ra,dec,w1mpro,w2mpro,w1sigmpro,cc_flags,ext_flg,nb
FROM allwise_p3as_psd
WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{CENTER},0.5))
AND w1mpro >= 7 AND w1mpro < 12 ORDER BY designation"""
# Additional faint neighbors cross-match around the field, not sample members.
NEIGHBORS = f"""SELECT designation,ra,dec,w1mpro,w2mpro,cc_flags,ext_flg,nb
FROM allwise_p3as_psd
WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{CENTER},0.51))
AND w1mpro < 16 ORDER BY designation"""
PHOT = WISE.replace("designation,ra,dec,w1mpro,w2mpro,w1sigmpro,cc_flags,ext_flg,nb", "designation,j_m_2mass,h_m_2mass,k_m_2mass,w1mpro,w2mpro,ph_qual")
BRIGHT = f"""SELECT designation,ra,dec,w1mpro,w2mpro
FROM allwise_p3as_psd
WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{CENTER},0.75))
AND w1mpro < 9 ORDER BY designation"""
GAIA = f"""SELECT source_id,ra,dec,pmra,pmdec,phot_g_mean_mag,teff_gspphot,parallax,ruwe
FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',{CENTER},0.51))
AND phot_g_mean_mag < 19 ORDER BY source_id"""

def now():
    return datetime.now(timezone.utc).isoformat()

def pinned(name, size, digest):
    with (HERE/'input_hashes.csv').open() as file:
        record=next(r for r in csv.DictReader(file) if r['field']==FIELD and r['name']==name)
    if (int(record['bytes']),record['sha256'])!=(size,digest):
        raise RuntimeError(f"{name}: archive response differs from frozen input")

def fetch_catalog(name, service, query):
    uri = service+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query,MAXREC=50000))
    start = now()
    with urlopen(uri,timeout=180) as r:
        data = r.read()
        if r.status!=200 or data.startswith((b'<',b'ERROR')):
            raise RuntimeError(data[:2000])
        hdr = r.headers
    lines = list(csv.DictReader(data.decode().splitlines()))
    if not lines or len(lines)>=50000: raise ValueError('empty or possibly truncated catalog')
    pinned(name,len(data),hashlib.sha256(data).hexdigest())
    (ROOT/'inputs'/f'{name}.csv').write_bytes(data)
    meta = dict(service=service,url=uri,query=query,retrieved_utc=start,rows=len(lines),maxrec=50000,
                http_status=200,content_type=hdr.get('Content-Type'),size=len(data),sha256=hashlib.sha256(data).hexdigest())
    (ROOT/'receipts'/f'{name}.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(name,len(lines),flush=True)

def fetch_cube():
    url = tile_url(TILE)
    start = now()
    with urlopen(Request(url,headers={'Range':'bytes=0-2879'}),timeout=90) as r:
        if r.status!=206 or not r.headers['Content-Range'].startswith('bytes 0-2879/'):
            raise IOError('first byte range mismatch')
        size = int(r.headers['Content-Range'].split('/')[1]); etag = r.headers['ETag']
        if len(r.read()) != 2880: raise IOError('incomplete initial block')
        modified = r.headers['Last-Modified']
    finalpath = ROOT/'localcache'/f'Npix{TILE}_cube.fits'
    with tempfile.NamedTemporaryFile(dir=finalpath.parent,prefix=f'Npix{TILE}.',suffix='.part',delete=False) as tmp:
        target=Path(tmp.name)
    step = 8*1024**2
    def download(lo):
        hi = min(lo+step,size)-1
        for attempt in range(4):
            try:
                with urlopen(Request(url,headers={'Range':f'bytes={lo}-{hi}','If-Match':etag}),timeout=100) as r:
                    if r.status!=206 or r.headers['Content-Range']!=f'bytes {lo}-{hi}/{size}' or r.headers['ETag']!=etag:
                        raise IOError('range, size, or ETag mismatch')
                    chunk=r.read()
                if len(chunk)!=hi-lo+1: raise IOError('truncated')
                with target.open('r+b') as f:
                    f.seek(lo); f.write(chunk)
                return dict(first=lo,last=hi,sha256=hashlib.sha256(chunk).hexdigest())
            except (OSError,ValueError):
                if attempt==3: raise
                time.sleep(2**attempt)
    with target.open('wb') as f: f.truncate(size)
    with ThreadPoolExecutor(10) as pool:
        ranges=list(pool.map(download,range(0,size,step)))
    with target.open('rb') as f: sha=hashlib.file_digest(f,'sha256').hexdigest()
    with target.open('rb') as f: md5=hashlib.file_digest(f,'md5').hexdigest()
    if '-' not in etag and md5!=etag.strip('"'): raise IOError('full cube MD5/ETag mismatch')
    pinned('cube',size,sha)
    target.replace(finalpath)
    meta = dict(url=url,tile=TILE,retrieved_utc=start,completed_utc=now(),size=size,etag=etag,
                last_modified=modified,sha256=sha,md5=md5,ranges=ranges)
    (ROOT/'receipts'/'cube.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('cube',sha,flush=True)

if __name__=='__main__':
    for folder in ['inputs','receipts','localcache']: (ROOT/folder).mkdir(exist_ok=True,parents=True)
    tasks = dict(cube=fetch_cube,allwise=lambda:fetch_catalog('allwise','https://irsa.ipac.caltech.edu/TAP/sync',WISE),
                 neighbors=lambda:fetch_catalog('neighbors','https://irsa.ipac.caltech.edu/TAP/sync',NEIGHBORS),
                 phot=lambda:fetch_catalog('phot','https://irsa.ipac.caltech.edu/TAP/sync',PHOT),
                 bright=lambda:fetch_catalog('bright','https://irsa.ipac.caltech.edu/TAP/sync',BRIGHT),
                 gaia=lambda:fetch_catalog('gaia','https://gea.esac.esa.int/tap-server/tap/sync',GAIA))
    with ThreadPoolExecutor(len(tasks)) as pool:
        jobs = {k:pool.submit(v) for k,v in tasks.items() if not sys.argv[1:] or k in sys.argv[1:]}
    for k,j in jobs.items():
        try: j.result()
        except Exception as err: print(k,repr(err),flush=True); raise
