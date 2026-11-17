"""Fetch FITS headers and the complete SCI, ERR and WHT planes of three F182M mosaics."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import time
from urllib.request import Request, urlopen

from astropy.io import fits


def read_range(url, start, end):
    for attempt in range(5):
        try:
            req=Request(f'{url}&range={start}-{end}',headers={'Range':f'bytes={start}-{end}'})
            with urlopen(req,timeout=95) as response:
                span=response.headers.get('Content-Range','')
                raw=response.read()
                if not span.startswith(f'bytes {start}-{end}/') or len(raw)!=end-start+1:
                    raise IOError('Incomplete response for '+req.full_url)
                return raw,int(span.split('/')[-1])
        except OSError:
            if attempt==4:raise
            time.sleep(2**attempt)


def product(visit, output):
    name=f'jw02768-o{visit:03d}_t001_nircam_clear-f182m_i2d.fits'
    url='https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/'+name
    tmp=output/(name+'.partial')
    with tmp.open('wb') as out:
        offset=0
        _,size=read_range(url,0,2879)
        while offset<size:
            header=bytearray()
            while True:
                block,_=read_range(url,offset+len(header),offset+len(header)+2879)
                header+=block
                if any(block[j:j+8]==b'END     ' for j in range(0,2880,80)):break
            h=fits.Header.fromstring(header.decode(),sep='')
            import math
            axes=h['NAXIS']
            elements=math.prod(h[f'NAXIS{i+1}'] for i in range(axes)) if axes else 0
            length=abs(h['BITPIX'])//8*(elements+h.get('PCOUNT',0))*h.get('GCOUNT',1)
            padded=(length+2879)//2880*2880
            out.seek(offset);out.write(header);offset+=len(header)
            if h.get('EXTNAME') in {'SCI','ERR','WHT'}:
                for start in range(offset,offset+padded,4*1024**2):
                    payload,_=read_range(url,start,min(start+4*1024**2,offset+padded)-1)
                    out.write(payload)
            offset+=padded
        out.truncate(size)
    tmp.replace(output/name)
    return describe(output/name)


def describe(path):
    extents=[]
    with fits.open(path) as hdus,path.open('rb') as inp:
        for h in hdus:
            info=h.fileinfo(); head=info['hdrLoc'];start=info['datLoc'];n=info['datSpan']
            inp.seek(head);hh=hashlib.sha256(inp.read(start-head)).hexdigest()
            record=dict(hdu=h.name,header_start=head,header_bytes=start-head,header_sha256=hh,
                        data_start=start,data_padded_bytes=n,retained=h.name in {'SCI','ERR','WHT'})
            if record['retained']:
                remaining=n;digest=hashlib.sha256()
                while remaining:
                    buf=inp.read(min(remaining,4*1024**2))
                    if not buf:raise EOFError(path)
                    digest.update(buf);remaining-=len(buf)
                record['data_padded_sha256']=digest.hexdigest()
            extents.append(record)
    return dict(filename=path.name,uri='mast:JWST/product/'+path.name,logical_bytes=path.stat().st_size,hdus=extents)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path)
    p.add_argument('--manifest',type=Path,required=True);p.add_argument('--existing',action='store_true')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    if a.existing:
        manifests=[describe(a.output/f'jw02768-o{v:03d}_t001_nircam_clear-f182m_i2d.fits') for v in (1,2,3)]
    else:
        with ThreadPoolExecutor(3) as ex:manifests=list(ex.map(lambda v:product(v,a.output),[1,2,3]))
    a.manifest.write_text(json.dumps(manifests,indent=2)+'\n')
