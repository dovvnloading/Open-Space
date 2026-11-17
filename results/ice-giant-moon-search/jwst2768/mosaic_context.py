"""Acquire original drizzle contribution masks alongside header-preserving sparse F182M mosaics."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path

from astropy.io import fits
from mosaic_fetch import read_range


def context(path, workers):
    url='https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/'+path.name
    with fits.open(path) as f:
        con=f['CON'];hdr=con.header.copy();item=con.fileinfo()
        start=item['datLoc'];pad=item['datSpan']
    shape=[hdr[f'NAXIS{i}'] for i in range(hdr['NAXIS'],0,-1)]
    if hdr['BITPIX']!=32 or hdr.get('BZERO',0)!=0 or hdr.get('BSCALE',1)!=1:
        raise ValueError('unexpected context dtype')
    from math import prod
    used=prod(shape)*4
    end=start+pad;stride=4*1024**2
    starts=range(start,end,stride)
    def block(offset):
        return read_range(url,offset,min(offset+stride,end)-1)[0]
    dest=path.with_suffix('.con.bin');tmp=path.with_suffix('.con.partial')
    digest=hashlib.sha256()
    with ThreadPoolExecutor(workers) as executor,tmp.open('wb') as out:
        for raw in executor.map(block,starts):out.write(raw);digest.update(raw)
    if tmp.stat().st_size!=pad:raise ValueError(tmp)
    tmp.replace(dest)
    return dict(source='mast:JWST/product/'+path.name,output=dest.name,hdu='CON',shape=shape,storage_dtype='>i4',
                data_start=start,payload_bytes=used,padded_bytes=pad,padded_sha256=digest.hexdigest())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);p.add_argument('manifest',type=Path)
    a=p.parse_args()
    with ThreadPoolExecutor(3) as ex:
        receipts=list(ex.map(lambda path:context(path,4),sorted(a.directory.glob('jw02768-o*_i2d.fits'))))
    if len(receipts)!=3:raise ValueError('expected three visits')
    a.manifest.write_text(json.dumps(receipts,indent=2)+'\n')
