"""Cache the official 2040 x 2040 detector wavelength/bandwidth maps used in DR1 and QR3."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits

S3 = 'https://nasa-irsa-spherex.s3.amazonaws.com/'
VERSION = 'cal-swcs-v5-2026-191'

def fetch(detector,cache,out):
    key=f'spectral_wcs/{VERSION}/{detector}/spectral_wcs_D{detector}_spx_{VERSION}.fits'
    path=cache/f'sci_D{detector}.fits'
    if not path.exists():
        with urlopen(S3+'dr1/'+key,timeout=120) as response:
            data=response.read()
        if len(data)!=33312960:
            raise ValueError('incomplete calibration image')
        path.write_bytes(data)
    data=path.read_bytes()
    records={}
    for release in ['dr1','qr3']:
        url=S3+release+'/'+key
        from urllib.request import Request
        with urlopen(Request(url,method='HEAD'),timeout=120) as response:
            records[release]=dict(url=url,etag=response.headers['ETag'].strip('"'))
    md5=hashlib.md5(data).hexdigest()
    if any(record['etag']!=md5 for record in records.values()):
        raise ValueError('release calibration mismatch')
    with fits.open(path) as h:
        print(detector,[(v.name,None if v.data is None else v.data.shape) for v in h],flush=True)
    records.update(sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))
    (out/f'science_wcs_D{detector}.json').write_text(json.dumps(records,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cache',type=Path)
    p.add_argument('out',type=Path)
    p.add_argument('detectors',type=int,nargs='+')
    a=p.parse_args()
    a.cache.mkdir(exist_ok=True,parents=True);a.out.mkdir(exist_ok=True,parents=True)
    with ThreadPoolExecutor(2) as pool:
        list(pool.map(lambda d:fetch(d,a.cache,a.out),a.detectors))
