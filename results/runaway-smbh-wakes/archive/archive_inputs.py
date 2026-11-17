"""Recreate the ancillary imaging inputs and acquisition receipt."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
from urllib import parse, request
import shutil
import time
from urllib.error import HTTPError

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

JADES = ('https://stpubdata.s3.amazonaws.com/mast/hlsp/jades/dr2/goods-s/images/'
         'hlsp_jades_jwst_nircam_goods-s-deep_f150w_v2.0_drz.fits')
FIELD = (53.16250302, -27.79141665)  # center of the original HST XDF mosaic
CONTROL = (40.4399,-8.3497)  # MAST HST GO-17301 target center

def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def range_read(url, start, end, etag):
    req = request.Request(url, headers={'Range':f'bytes={start}-{end}', 'If-Match':etag})
    for attempt in range(5):
        try:
            result=request.urlopen(req, timeout=120)
            break
        except HTTPError as exc:
            if attempt==4 or exc.code not in (429,502,503,504):raise
            time.sleep(3*(attempt+1))
    with result as r:
        if r.status != 206 or not r.headers['Content-Range'].startswith(f'bytes {start}-{end}/'):
            raise RuntimeError(f'Invalid HTTP range {r.headers.get("Content-Range")}')
        data = r.read()
        if len(data) != end-start+1:
            raise RuntimeError('Incomplete range')
    return data

def jades(out, center=FIELD, nside=1024, prefix="jades_f150w"):
    with request.urlopen(request.Request(JADES, method='HEAD'),timeout=60) as r:
        etag = r.headers['ETag']
        filesize = int(r.headers['Content-Length'])
    offset = 0
    receipt = {'url':JADES,'etag':etag,'object_bytes':filesize,'planes':[]}
    while offset < filesize:
        header_data = range_read(JADES, offset, min(offset+57599,filesize-1), etag)
        end = next(i for i in range(0,len(header_data),80) if header_data[i:i+8]==b'END     ')
        skip = math.ceil((end+80)/2880)*2880
        h = fits.Header.fromstring(header_data[:end+80].decode('ascii'), sep='')
        offset += skip
        plane = h.get('EXTNAME','PRIMARY').strip()
        if h['NAXIS'] == 0:
            continue
        nx,ny = h['NAXIS1'],h['NAXIS2']
        bpp = abs(h['BITPIX']) // 8
        if plane in ('SCI','WHT','ERR'):
            wcs = WCS(h).celestial
            x,y = wcs.all_world2pix([center],0)[0]
            # Rounded fixed physical center. Crop science and registered error on identical pixels.
            if plane=='SCI':
                x0,y0 = round(x)-nside//2,round(y)-nside//2
            start = offset + y0*nx*bpp
            count = nside*nx*bpp
            chunk = nx*bpp*32
            def read(i):
                a,b = start+i,min(start+i+chunk,start+count)-1
                return range_read(JADES,a,b,etag)
            with ThreadPoolExecutor(max_workers=5) as ex:
                data = b''.join(ex.map(read,range(0,count,chunk)))
            a = np.frombuffer(data,dtype='>f4').reshape(nside,nx)[:,x0:x0+nside].copy()
            cutwcs = wcs.slice((slice(y0,y0+nside),slice(x0,x0+nside)))
            ch = cutwcs.to_header()
            ch['BUNIT']=h.get('BUNIT',''); ch['FILTER']='F150W'
            dst = out/f'{prefix}_{plane.lower()}.fits'
            fits.writeto(dst,a,ch,overwrite=True)
            receipt['planes'].append(dict(hdu=plane,xywh=[int(x0),int(y0),nside,nside],header_offset=offset-skip,
                                    first_byte=start,byte_count=count,row_stride=nx*bpp,range_block=chunk,
                                    rows_sha256=hashlib.sha256(data).hexdigest(),file=dst.name,sha256=digest(dst)))
        offset += math.ceil(ny*nx*bpp/2880)*2880
        if plane=='WHT':
            break
    return receipt

def other(out, name, url, existing):
    dest=out/name
    if not (existing and dest.exists()):
        for attempt in range(4):
            try:
                partial = dest.with_suffix('.part')
                with request.urlopen(url, timeout=100) as rx, partial.open('wb') as fd:
                    shutil.copyfileobj(rx,fd)
                partial.replace(dest)
                break
            except HTTPError as exc:
                if exc.code not in (429,502,503,504) or attempt==3:
                    raise
                time.sleep(5*(attempt+1))
    with fits.open(dest) as hdus:
        shape=[getattr(h,'shape',None) for h in hdus]
    return dict(file=name, url=url, sha256=digest(dest),size=dest.stat().st_size,hdu_shapes=shape)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    p.add_argument('--existing',action='store_true',help='reuse local full files, recomputing hashes')
    p.add_argument('--receipt',type=Path,required=True)
    args=p.parse_args()
    out=args.directory;out.mkdir(parents=True,exist_ok=True)
    catalog='hlsp_uvudf_hst_wfc3-uvis_udf-epoch3_multi_v2.0_cat.fits'
    urls={catalog:'https://archive.stsci.edu/missions/hlsp/uvudf/v2.0/'+catalog}
    for band in ('f606w','f814w'):
        for typ in ('sci','wht'):
            fn=f'hlsp_xdf_hst_acswfc-60mas_hudf_{band}_v1_{typ}.fits'
            urls[fn]='https://archive.stsci.edu/hlsps/xdf/'+fn
    for obs in ['jety02010','jety02020']:
        fn=obs+'_drc.fits'
        urls[fn]='https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/'+fn
    for field,center in [('xdf',FIELD),('rbh1',CONTROL)]:
        for release in ['ls-dr10','ls-dr11']:
            fn=f'legacy_{field}_{release}.fits'
            urls[fn]='https://www.legacysurvey.org/viewer/fits-cutout?'+parse.urlencode(
                dict(ra=center[0],dec=center[1],layer=release,pixscale=.262,size=256,bands='grz'))
    def fetch(item):
        fn,url=item
        try: return other(out,fn,url,args.existing)
        except Exception as exc: return dict(file=fn,url=url,error=str(exc))
    with ThreadPoolExecutor(max_workers=4) as ex: rec=list(ex.map(fetch,urls.items()))
    rec.append(jades(out))
    args.receipt.write_text(json.dumps(rec,indent=2)+'\n')

if __name__=='__main__':
    main()
