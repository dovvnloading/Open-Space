"""Fetch original VLASS sky pixels for WISE J0458+6434 using public CADC FITS byte ranges."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from replay import predicted


def ranged(url, lo, hi):
    headers = {"Range": f"bytes={lo}-{hi}"}
    with urlopen(Request(url,headers=headers), timeout=120) as response:
        if response.status != 206 or not response.headers["Content-Range"].startswith(f"bytes {lo}-{hi}/"):
            raise IOError("Server returned a different range")
        result = response.read()
        if len(result) != hi-lo+1:
            raise IOError("Truncated range")
        return result


def extract(epoch, out, side):
    uri = f"VLASS{epoch}.1.ql.T26t04.J045918+643000.10.2048.v1.I.iter1.image.pbcor.tt0.subim.fits"
    url = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/VLASS/" + quote(uri,safe="")
    first = ranged(url,0,115199)
    end = next(i+80 for i in range(0,len(first),80) if first[i:i+8] == b"END     ")
    header_size = ((end+2879)//2880)*2880
    hdr = fits.Header.fromstring(first[:end].decode("ascii"))
    if hdr["BITPIX"] != -32 or hdr["NAXIS3"] != 1 or hdr["NAXIS4"] != 1:
        raise ValueError("Unexpected dimensions or pixel type")
    ra, dec = predicted(hdr["DATE-OBS"])
    cx, cy = np.round(WCS(hdr).celestial.world_to_pixel_values(ra,dec)).astype(int)
    x0, y0 = cx-side//2, cy-side//2
    width = hdr["NAXIS1"]
    lo = header_size + y0*width*4
    hi = lo+side*width*4-1
    payload = ranged(url,lo,hi)
    data = np.frombuffer(payload,dtype=">f4").reshape(side,width)[:,x0:x0+side]
    hdr["CRPIX1"] -= x0
    hdr["CRPIX2"] -= y0
    fn=out/f"VLASS{epoch}.fits"
    fits.PrimaryHDU(data=data[None,None,:,:],header=hdr).writeto(fn,overwrite=True)
    return dict(epoch=epoch,url=url,side=side,header_block_bytes=header_size,full_row_byte_range=[int(lo),int(hi)],
                pixel_sha256=hashlib.sha256(np.asarray(data,dtype=">f4").tobytes()).hexdigest(),
                original_rows_sha256=hashlib.sha256(payload).hexdigest())


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("output",type=Path)
    p.add_argument("--side",type=int,default=91)
    p.add_argument("epochs",nargs="*",type=int,default=[1,2,3,4])
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    receipts=[extract(i,args.output,args.side) for i in args.epochs]
    (args.output/"download.json").write_text(json.dumps(receipts,indent=2)+"\n")
