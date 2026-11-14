"""Augment 97-pixel sky stamps with local exposure map and local hashes."""
from swift_stamps import D,BASE,POS,RECEIPTS
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from io import BytesIO
from urllib.request import urlopen
import json,hashlib,gzip
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
def add(p):
 o=json.loads(p.read_text());q=p.with_suffix('.exp.npz');meta=p.with_suffix('.expmeta')
 if not o['planes'] or q.exists():return
 with np.load(p.with_suffix('.npz')) as imgs:
    if not any(z['counts']>0 for z in o['planes']):return
 url=o['url'].replace('_sk.img','_ex.img')
 try:
    expected=RECEIPTS[o['obsid']]['sha256_exposure']
    if not expected:return
    with urlopen(url,timeout=35) as response:
        b=response.read()
    if hashlib.sha256(b).hexdigest()!=expected:raise ValueError('Archive exposure checksum mismatch')
    a={}
    with fits.open(BytesIO(gzip.decompress(b))) as hd:
        for z in o['planes']:
            i=z['idx'];h=hd[i]
            w=WCS(h.header).celestial
            xy=w.all_world2pix([POS],0)[0]
            if np.max(np.abs(xy-z['target']))>.1 or list(h.data.shape)!=list(reversed(z['dims'])):continue
            ix,iy=z['ci']
            a['exp'+str(i)]=h.data[iy-48:iy+49,ix-48:ix+49].copy()
    if a:np.savez_compressed(q,**a);meta.write_text(hashlib.sha256(b).hexdigest())
 except Exception as e:raise RuntimeError(f'Could not recover Swift exposure {p.name}') from e
if __name__=='__main__':
 with ThreadPoolExecutor(12) as ex:list(ex.map(add,(D/'swift_stamps').glob('*.json')))
