"""Pull small stamps from authenticated-by-hash public QDR image products."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
import gzip
import json
from pathlib import Path
from urllib.request import Request, urlopen
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from extract import BASE, HashingReader

R=27

def extract_stamps(item):
    ep,band,kind,pos,cache,out=item
    filename=f'nexus_central_deep_ep{ep:02d}_{band}'+('_60mas' if band=='f444w' else '')+f'_i2d_{kind}.fits.gz'
    url=BASE+filename
    dst=out/f'ep{ep:02d}_{band}_{kind}.npz'
    if dst.exists():
        old=np.load(dst)
        if np.array_equal(old['id'],pos.id.to_numpy(dtype=str)): return
    src=cache/filename
    with closing(src.open('rb') if src.exists() else urlopen(Request(url,headers={'Range':'bytes=0-'}),timeout=180)) as raw:
        wrapped=HashingReader(raw)
        with gzip.GzipFile(fileobj=wrapped) as im:
            h=fits.Header.fromfile(im)
            nx,ny=h['NAXIS1'],h['NAXIS2']
            assert h['BITPIX']==-32
            w=WCS(h)
            px,py=w.all_world2pix(pos.ra_deg,pos.dec_deg,0)
            xi,yi=np.rint(px).astype(int),np.rint(py).astype(int)
            a=np.full((len(pos),2*R+1,2*R+1),np.nan,'f4')
            # Read only stripes that intersect a stamp, traversing gzip in order.
            start=im.tell()
            for y0 in sorted({(yy//128)*128 for yy0 in yi for yy in [yy0-R,yy0,yy0+R] if 0<=yy<ny}):
                im.seek(start+4*nx*y0)
                y1=min(ny,y0+128)
                data=np.frombuffer(im.read((y1-y0)*nx*4),'>f4').reshape(-1,nx)
                for j in np.flatnonzero((yi+R>=y0)&(yi-R<y1)&(xi-R>=0)&(xi+R<nx)):
                    lo=max(yi[j]-R,y0); hi=min(yi[j]+R+1,y1)
                    a[j,lo-yi[j]+R:hi-yi[j]+R]=data[lo-y0:hi-y0,xi[j]-R:xi[j]+R+1]
            while im.read(4*1024*1024): pass
        while wrapped.read(4*1024*1024): pass
        if not src.exists() and raw.headers.get('Content-Length'):
            assert wrapped.size==int(raw.headers['Content-Length'])
    np.savez_compressed(dst,stamps=a,id=pos.id.to_numpy(dtype=str),fx=px-xi,fy=py-yi,
                        njy_per_unit=h['PIXAR_SR']*1e15)
    dst.with_suffix('.json').write_text(json.dumps(dict(epoch=ep,filter=band,product=kind,
        original_url=url,compressed_bytes=wrapped.size,compressed_sha256=wrapped.sha.hexdigest()),indent=2)+'\n')
    print(dst.name,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('positions',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--cache',type=Path,default=Path('.'));p.add_argument('--epochs',nargs='+',type=int,default=list(range(1,9)))
    p.add_argument('--jobs',type=int,default=4);a=p.parse_args()
    pos=pd.read_csv(a.positions);a.out.mkdir(exist_ok=True,parents=True)
    jobs=[(e,b,k,pos,a.cache,a.out) for e in a.epochs for b in ['f444w','f200w'] for k in ['data','error']]
    with ThreadPoolExecutor(a.jobs) as ex: list(ex.map(extract_stamps,jobs))
