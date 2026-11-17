"""Find native Q1 VIS individual-exposure pixels for a specified sky segment.

IRSA IBE returns primary-image cutouts only; the SWL-DET exposes detector
quadrants in an uncompressed FITS HDU sequence. Check all HDU names, shapes,
stride and returned byte ranges before interpreting the selected stripe.
"""
import argparse
import csv
import hashlib
import json
import math
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

import numpy as np
from astropy.io.fits import Header
from astropy.wcs import WCS

from acquire import S3, sha
from filament import measure

BLOCK=2880


def get_range(url, start, end, logs, cache):
    key=hashlib.sha256((url+f'#{start}:{end}').encode()).hexdigest()
    dest=cache/key
    if dest.exists():
        raw=dest.read_bytes()
    else:
        full=url+f'?range_start={start}&range_end={end}'
        with urlopen(Request(full,headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'}),timeout=45) as f:
            if f.status!=206 or f.headers['Content-Range'].split('/')[0]!=f'bytes {start}-{end}':
                raise ValueError('unexpected byte range')
            raw=f.read()
        if len(raw)!=end-start+1:raise ValueError('incomplete range')
        dest.write_bytes(raw)
    logs.append(dict(start=start,end=end,sha256=sha(raw)))
    return raw


def header(raw):
    for i in range(0,len(raw),80):
        if raw[i:i+8]==b'END     ':
            n=((i+80+BLOCK-1)//BLOCK)*BLOCK
            return Header.fromstring(raw[:i+80].decode(),sep=''),n
    raise ValueError('missing END')


def sky_files(ra,dec,cache):
    url='https://irsa.ipac.caltech.edu/SIA?'+urlencode(dict(POS=f'CIRCLE {ra:.9f} {dec:.9f} 0.0005',
        COLLECTION='euclid_DpdVisCalibratedQuadFrame',MAXREC=1500))
    path=cache/(sha(url.encode())+'.xml')
    if not path.exists():
        with urlopen(url,timeout=60) as f: raw=f.read()
        root=ET.fromstring(raw)
        if not all(i.get('value')=='OK' for i in root.findall('.//{*}INFO[@name="QUERY_STATUS"]')):
            raise ValueError('SIA error or overflow')
        path.write_bytes(raw)
    raw=path.read_bytes()
    e=ET.fromstring(raw);heads=[x.get('name') for x in e.findall('.//{*}FIELD')]
    rows=[dict(zip(heads,[x.text for x in t.findall('{*}TD')])) for t in e.findall('.//{*}TR')]
    keys=sorted({json.loads(r['cloud_access'])['aws']['key'] for r in rows if 'SWL-DET' in r['access_url']})
    return keys,dict(url=url,sha256=sha(raw),bytes=len(raw),rows=len(rows),science_frames=len(keys))


def measure_frame(key,ends,cache):
    url=S3+key;logs=[];rec=dict(key=key,ranges=logs);center=ends.mean(axis=0)
    try:
        raw=get_range(url,0,57599,logs,cache)
        primary,primary_size=header(raw);first,first_size=header(raw[primary_size:])
        nx,ny=first['NAXIS1'],first['NAXIS2'];stride=first_size+math.ceil(nx*ny*4/BLOCK)*BLOCK
        if first['EXTNAME']!='1-1.E.SCI' or first['BITPIX']!=-32 or (nx,ny)!=(2048,2066):
            raise ValueError('nonstandard first quadrant')
        memo={0:(first,primary_size+first_size)}
        def get(q):
            if q not in memo:
                start=primary_size+q*stride
                raw=get_range(url,start,start+first_size-1,logs,cache)
                h,n=header(raw)
                expected=f'{q//72+1}-{q//12%6+1}.{"EFGH"[q//3%4]}.SCI'
                if not (h['EXTNAME']==expected and h['NAXIS1']==nx and h['NAXIS2']==ny and h['BITPIX']==-32 and n==first_size):
                    raise ValueError('unexpected quadrant signature '+str(q))
                memo[q]=(h,start+n)
            return memo[q]
        def position(q):
            h,_=get(q)
            return WCS(h).all_pix2world([[nx/2,ny/2]],0)[0]
        anchors=[0,3,6,9,12,72]
        pos={q:position(q) for q in anchors}
        dra=np.cos(np.deg2rad(center[1]))
        def dist(q):
            p=pos[q%12]+(q//72)*(pos[72]-pos[0])+(q//12%6)*(pos[12]-pos[0])
            return np.linalg.norm((p-center)*[dra,1.])
        indices=sorted(range(0,432,3),key=dist)
        # Six approximate nearest quadrants suffice for the regular camera grid.
        for idx in indices[:6]:
            h,offset=get(idx);w=WCS(h)
            try:
                xy=w.all_world2pix(ends,0)
            except Exception:
                continue
            margin=65
            if np.isnan(xy).any() or xy[:,0].min()<=margin or xy[:,1].min()<=margin or xy[:,0].max()>=nx-margin or xy[:,1].max()>=ny-margin:
                continue
            a=int(np.floor(xy[:,1].min()))-margin;b=int(np.ceil(xy[:,1].max()))+margin
            data=get_range(url,offset+a*nx*4,offset+(b+1)*nx*4-1,logs,cache)
            stripe=np.frombuffer(data,dtype='>f4').reshape(b-a+1,nx)
            scale=3600*np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))
            z=measure(stripe,xy[0]-[0,a],xy[1]-[0,a],scale,mask=np.isfinite(stripe))
            rec.update(status='measured', quadrant=h['EXTNAME'],xy=xy.tolist(),y_first=a,y_last=b,scale=scale,**z)
            return rec
        rec['status']='no_usable_interior_quadrant'
    except (OSError,ValueError,KeyError) as exc:
        rec.update(status='error',error=repr(exc))
    return rec


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--candidates',type=Path,required=True)
    ap.add_argument('--id',required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--jobs',type=int,default=4);a=ap.parse_args()
    a.cache.mkdir(parents=True,exist_ok=True);a.out.mkdir(parents=True,exist_ok=True)
    row=next(r for r in csv.DictReader(a.candidates.open()) if r['id']==a.id)
    ends=np.array([[float(row[f'ra{i}_deg']),float(row[f'dec{i}_deg'])] for i in (0,1)])
    keys,sia=sky_files(*ends.mean(axis=0),a.cache)
    with ThreadPoolExecutor(a.jobs) as pool:results=list(pool.map(lambda k:measure_frame(k,ends,a.cache),keys))
    (a.out/(a.id+'.json')).write_text(json.dumps(dict(id=a.id,sia=sia,frames=results),indent=1)+'\n')
    print(a.id,[(r['status'],r.get('central_ridge_contrast_sigma')) for r in results],flush=True)


if __name__=='__main__':main()
