"""Recover parents lacking usable Zenodo cutouts directly from the IRSA Q1 MER mosaics."""
import argparse
import csv
import json
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

from followup import cut,ring_profile,render
from photometry import keys_by_band,read,read_primary
from search import HERE,digest,find_rings,world

def recover(row,tile):
    keys=keys_by_band(tile)
    arrays=[];receipt=[]
    for name in ('image','rms'):
        key=keys['VIS',name]
        url='https://irsa.ipac.caltech.edu/ibe/data/euclid/'+key+f"?center={float(row['right_ascension']):.8f},{float(row['declination']):.8f}deg&size=32arcsec&gzip=false"
        raw=read(url)
        header,a=read_primary(raw);arrays.append(a)
        receipt.append(dict(product=name,url=url,sha256=digest(raw),bytes=len(raw)))
    flux,rms=arrays
    ny,nx=flux.shape[0]//2,flux.shape[1]//2
    a=flux[:2*ny,:2*nx].astype('f4');rms=rms[:2*ny,:2*nx]
    mask=np.isfinite(a)&(a!=0)&np.isfinite(rms)&(rms>0)&(rms<1e8)
    a[~mask]=0
    a=a.reshape(ny,2,nx,2).mean(axis=(1,3))
    mask=mask.reshape(ny,2,nx,2).all(axis=(1,3))
    peaks,n=find_rings(a,mask)
    positions=[]
    for j,(_,x,y,r,_,_) in enumerate(peaks):
        m=ring_profile(cut(a,x,y),r)
        ra,dec=world(header,2*x+.5,2*y+.5)
        positions.append(dict(ra=ra,dec=dec,r=r*.2,**m))
        if m['pass_vis']:
            render(cut(a,x,y),HERE/'stamps'/f"{row['id_str']}_fallback_{j}.png")
    return dict(id_str=row['id_str'],tile=tile,status='examined' if n else 'insufficient_coverage',
                eligible=n,triggers=len(peaks),pass_vis=sum(p['pass_vis'] for p in positions),positions=positions,inputs=receipt)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('original_catalog',type=Path)
    args=p.parse_args()
    tile_by_id={r['id_str']:r['tile_index'] for r in csv.DictReader(args.original_catalog.open())}
    catalogue={r['id_str']:r for r in csv.DictReader((HERE.parent/'catalog/slde_screen.csv').open())}
    missing=[r['id_str'] for r in json.loads((HERE/'unmodelled.json').read_text()) if r['status']!='examined']
    def run(i):
        try:
            return recover(catalogue[i],tile_by_id[i])
        except (OSError,ValueError,KeyError,AssertionError) as e:
            return dict(id_str=i,status='error',error=repr(e))
    with ThreadPoolExecutor(5) as pool:result=list(pool.map(run,missing))
    (HERE/'fallback.json').write_text(json.dumps(result,indent=2)+'\n')
    print([(r['id_str'],r['status'],r.get('pass_vis')) for r in result])
if __name__=='__main__':main()
