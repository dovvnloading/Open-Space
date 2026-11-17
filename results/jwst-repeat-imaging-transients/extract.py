"""Stream NEXUS QDR primary HDUs into WCS-preserving cutouts and hash inputs."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from pathlib import Path
from urllib.request import Request, urlopen
import argparse
import csv
import gzip
import hashlib
import json
import os
import numpy as np

ROOT=Path(__file__).resolve().parent
BASE='https://ariel.astro.illinois.edu/nexus/qdr/nircam/'
BOUNDS=(7680,6976,9728,9024)
class HashingReader:
    def __init__(self,stream):
        self.stream=stream
        self.sha=hashlib.sha256()
        self.size=0
    def read(self,n=-1):
        block=self.stream.read(n)
        self.sha.update(block)
        self.size+=len(block)
        return block

def read_header(stream):
    cards=[]
    while True:
        block=stream.read(2880)
        if len(block)!=2880: raise ValueError("Truncated FITS header")
        cards.extend(block[i:i+80] for i in range(0,2880,80))
        if any(c.startswith(b"END     ") for c in cards): break
    values={c[:8].decode().strip():c[10:].split(b'/')[0].decode().strip() for c in cards if c[8:10]==b'= '}
    return cards,values

def extract(ep,band,kind, cache, out, bounds=BOUNDS):
    filename=f'nexus_central_deep_ep{ep:02d}_{band}'+('_60mas' if band=='f444w' else '')+f'_i2d_{kind}.fits.gz'
    remote=BASE+filename
    cached=cache/filename
    url=remote
    dest=out/f'ep{ep:02d}_{band}_{kind}.fits'
    out.mkdir(parents=True,exist_ok=True)
    with closing(cached.open('rb') if cached.exists() else urlopen(Request(url,headers={'Range':'bytes=0-'}),timeout=120)) as raw:
        total=os.stat(cached).st_size if cached.exists() else int(raw.headers['Content-Length'])
        source=HashingReader(raw)
        with gzip.GzipFile(fileobj=source) as f:
            cards,h=read_header(f)
            nx,ny=int(h['NAXIS1']),int(h['NAXIS2'])
            if h['BITPIX']!='-32' or h['NAXIS']!='2': raise ValueError(h)
            assert abs(float(h['CRVAL1'])-268.4625)<1e-8 and abs(float(h['CRVAL2'])-65.199166666667)<1e-8
            assert float(h['PC1_1'])==-1 and all(float(h.get(k,0))==0 for k in ('PC1_2','PC2_1')) and float(h.get('PC2_2',1))==1
            ratio=(0.06/3600)/float(h['CDELT1']); assert abs(ratio-round(ratio))<1e-8 and abs(float(h['CDELT1'])-float(h['CDELT2']))<1e-12
            x0,y0,x1,y1=bounds
            def index(v,crpix,reference):
                edge=(v+0.5-reference)*ratio+crpix-0.5
                if abs(round(edge)-edge)>1e-5: raise ValueError(('pixel edge',edge))
                return round(edge)
            ax0,ax1=(index(v,float(h['CRPIX1']),8682.5) for v in (x0,x1))
            ay0,ay1=(index(v,float(h['CRPIX2']),8028.5) for v in (y0,y1))
            updates={'CRPIX1':float(h['CRPIX1'])-ax0,'CRPIX2':float(h['CRPIX2'])-ay0,'NAXIS1':ax1-ax0,'NAXIS2':ay1-ay0}
            with dest.with_suffix('.tmp').open('wb') as cut:
                for c in cards:
                    key=c[:8].decode().strip()
                    if key in updates: c=(f'{key:8}= {updates[key]:>20}').ljust(80).encode()
                    cut.write(c)
                remaining=ay0*nx*4
                while remaining:
                    b=f.read(min(4*1024*1024,remaining))
                    if not b: raise ValueError('unexpected end')
                    remaining-=len(b)
                for j in range(ay0,ay1,32):
                    rows=min(32,ay1-j)
                    block=f.read(nx*rows*4)
                    a=np.frombuffer(block,dtype='>f4').reshape(rows,nx)
                    cut.write(a[:,ax0:ax1].tobytes())
                cut.write(bytes((-cut.tell())%2880))
        while source.read(4*1024*1024): pass
        if total!=source.size: raise ValueError((total,source.size))
    dest.with_suffix('.tmp').replace(dest)
    rec={'epoch':ep,'filter':band,'product':kind,'original_url':url,'compressed_bytes':total,'compressed_sha256':source.sha.hexdigest(),'parent_x0':ax0,'parent_x1':ax1,'parent_y0':ay0,'parent_y1':ay1,'cutout_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'cutout_file':dest.name}
    dest.with_suffix('.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(dest.name,source.size,source.sha.hexdigest(),flush=True)
    return rec

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--epochs',nargs='+',type=int,default=[7,8])
    parser.add_argument('--bands',nargs='+',choices=['f444w','f200w'],default=['f444w','f200w'])
    parser.add_argument('--cache',type=Path,default=ROOT)
    parser.add_argument('--out',type=Path,default=ROOT/'cuts')
    parser.add_argument('--jobs',type=int,default=4)
    parser.add_argument('--bounds',type=int,nargs=4,default=BOUNDS,metavar=('X0','Y0','X1','Y1'))
    args=parser.parse_args()
    with ThreadPoolExecutor(args.jobs) as executor:
        jobs=[executor.submit(extract,i,j,k,args.cache,args.out,args.bounds) for j in args.bands for i in args.epochs for k in ['data','error']]
        for job in jobs: job.result()
