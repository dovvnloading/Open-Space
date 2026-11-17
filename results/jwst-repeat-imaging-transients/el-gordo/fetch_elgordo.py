"""Read rectangular ranges of original STScI public i2d mosaics with byte receipts."""
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import sys
import time
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

HERE = Path(__file__).resolve().parent
CACHE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("cluster-cache")
CACHE.mkdir(exist_ok=True, parents=True)
CENTER = (15.735, -49.252)
SIZE_ARCSEC = 76
OBS = ["jw01176-o241_t012", "jw06882-o051_t051"]
BANDS = ["f444w", "f150w"]


def fetch(args):
    obs, band = args
    root = f"{obs}_nircam_clear-{band}"
    program, visit = obs.split("_")[0].split("-")
    url = f"https://stpubdata.s3.amazonaws.com/jwst/public/{program}/L3/t/{visit}/{root}_i2d.fits"
    parts = []
    def get(start, size):
        expected = f"bytes {start}-{start + size - 1}/"
        for attempt in range(4):
            try:
                request = Request(url, headers={"Range": f"bytes={start}-{start + size - 1}"})
                with urlopen(request, timeout=50) as response:
                    assert response.status == 206 and response.headers["Content-Range"].startswith(expected)
                    etag = response.headers['ETag']
                    body = response.read()
                assert len(body) == size
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(1+attempt)
        parts.append(dict(first=start, length=size, sha256=sha256(body).hexdigest(), etag=etag))
        return body

    offset = 0
    out = []
    window = None
    for index in range(15):
        blocks = bytearray()
        while True:
            block = get(offset + len(blocks), 2880 * 4)
            blocks.extend(block)
            cards = [blocks[i:i+80] for i in range(0, len(blocks), 80)]
            ends = [i for i, card in enumerate(cards) if card[:8] == b'END     ']
            if ends:
                length = 2880 * ((ends[0] * 80 + 80 + 2879)//2880)
                header = fits.Header.fromstring(blocks[:length].decode("ascii"), sep="")
                break
        pos = offset + length
        dims = tuple(header[f"NAXIS{i+1}"] for i in range(header["NAXIS"]))
        count = int(np.prod(dims, dtype=np.int64)) if dims else 0
        element_bytes = abs(header["BITPIX"])//8
        offset = pos + 2880 * ((count * element_bytes + 2879)//2880)
        name = header.get("EXTNAME", "PRIMARY")
        if index == 0:
            out.append(fits.PrimaryHDU(header=header))
        if name == "SCI":
            wcs = WCS(header, naxis=2)
            cx,cy = wcs.world_to_pixel_values(*CENTER)
            half = int(np.ceil(SIZE_ARCSEC/(3600*min(proj_plane_pixel_scales(wcs)))*np.sqrt(2)*0.5)) + 20
            x, y = int(round(float(cx))), int(round(float(cy)))
            window = max(0,x-half),min(dims[0],x+half+1),max(0,y-half),min(dims[1],y+half+1)
        if name in ('SCI','ERR','WHT','CON'):
            x0,x1,y0,y1 = window
            dtype={-32:'>f4',-64:'>f8',32:'>i4',16:'>i2',8:'u1'}[header["BITPIX"]]
            planes=[]
            for j in range(dims[2] if len(dims)==3 else 1):
                start = pos + (j*dims[1]+y0)*dims[0]*element_bytes
                n = (y1-y0)*dims[0]*element_bytes
                # Small requests support retries without retaining a full mosaic in memory.
                chunk=8*1024**2
                def stripe(i): return get(start+i,min(chunk,n-i))
                with ThreadPoolExecutor(4) as workers:
                    b=b''.join(workers.map(stripe,range(0,n,chunk)))
                arr=np.frombuffer(b,dtype=dtype).reshape(y1-y0,dims[0])[:,x0:x1].copy()
                planes.append(arr)
            sliced = header.copy()
            for key,value in [('CRPIX1',x0),('CRPIX2',y0)]:
                if key in sliced: sliced[key]-=value
            out.append(fits.ImageHDU(planes[0] if len(dims)==2 else np.stack(planes),header=sliced))
        if set(h.name for h in out) >= {'PRIMARY','SCI','ERR','CON','WHT'}: break
    assert set(h.name for h in out) >= {'PRIMARY','SCI','ERR','CON','WHT'}
    fits.HDUList(out).writeto(CACHE/(root+'.fits'),overwrite=True,output_verify="silentfix")
    print(root,window,out[1].data.shape,flush=True)
    return dict(root=root,url=url,window_xyxy=[x0,y0,x1,y1],ranges=sorted(parts,key=lambda p:p['first']),cutout_sha256=sha256((CACHE/(root+'.fits')).read_bytes()).hexdigest())


if __name__ == '__main__':
    entries=list(ThreadPoolExecutor(4).map(fetch,[(o,b) for o in OBS for b in BANDS]))
    (HERE/'input_ranges.json').write_text(json.dumps(entries,indent=1)+'\n')
