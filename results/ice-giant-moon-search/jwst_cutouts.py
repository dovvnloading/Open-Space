"""Download native CAL SCI, ERR and DQ windows from JWST public S3 FITS files.

The cutout preserves SIP coefficients and shifts CRPIX into the cropped pixel frame.
The provenance ledger hashes source bytes by original FITS file offset, and the
result ledger hashes each output. Only chips intersecting the square of interest
are downloaded. Disjoint chips are listed in cutouts.csv with zero bytes.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import io
from pathlib import Path
import time
from urllib.request import Request, urlopen
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def make_cut(row, half, dest):
    filename=row["productFilename"]
    url=("https://stpubdata.s3.amazonaws.com/jwst/public/jw"+filename[2:7]+"/"+
         filename.split("_")[0]+"/"+filename)
    receipts=[]
    def get(lo,n):
        hi=lo+n-1
        for a in range(4):
            try:
                with urlopen(Request(url+f"?bytes={lo}-{hi}",headers={"Range":f"bytes={lo}-{hi}"}),timeout=120) as r:
                    if r.status !=206 or r.headers.get("Content-Range")!=f"bytes {lo}-{hi}/{row['size']}":
                        raise OSError("invalid content-range")
                    b=r.read(n+1)
                if len(b)!=n: raise OSError("short response")
                receipts.append((filename,lo,n,hashlib.sha256(b).hexdigest()))
                return b
            except OSError:
                if a==3:raise
                time.sleep(2**a)
    raw=io.BytesIO(get(0,65536))
    ph=fits.Header.fromfile(raw)
    h=fits.Header.fromfile(raw)
    mid=ph["EXPMID"];sx,sy=h["NAXIS1"],h["NAXIS2"]
    x,y=WCS(h,naxis=2).all_world2pix(h["MT_RA"],h["MT_DEC"],0)
    xmin,xmax=max(0,int(x)-half),min(sx,int(x)+half+1)
    ymin,ymax=max(0,int(y)-half),min(sy,int(y)+half+1)
    base=(filename, ph["FILTER"],ph["VISIT_ID"],mid,round(float(x),3),round(float(y),3))
    if xmax<=xmin or ymax<=ymin:return (*base,"","","",0,""),receipts
    ph["CUTX0"]=xmin;ph["CUTY0"]=ymin
    ph["CUTSRC"]=filename
    ph.add_history("SCI, ERR and DQ windows from native calibrated file; CRPIX transformed to cutout.")
    hdus=[fits.PrimaryHDU(header=ph)];p=raw.tell()
    for i,kind in enumerate(["SCI","ERR","DQ"]):
        if i:
            h=fits.Header.fromfile(io.BytesIO(get(p,2880)));p+=2880
        if h["EXTNAME"]!=kind or abs(h['BITPIX'])!=32:raise ValueError("unexpected FITS HDU")
        b=get(p+ymin*sx*4,(ymax-ymin)*sx*4)
        arr=np.frombuffer(b,">f4" if h["BITPIX"]<0 else ">i4").reshape(-1,sx)[:,xmin:xmax]
        p+=((h.data_size+2879)//2880)*2880
        if kind=="DQ":
            if h.get("BZERO",0)==2**31:
                arr=np.bitwise_xor(arr.astype(np.uint32),np.uint32(2**31))
                del h["BZERO"]
                h.pop('BSCALE',None)
            else:arr=arr.astype(np.int32)
        if "CRPIX1" in h:h["CRPIX1"]-=xmin
        if "CRPIX2" in h:h["CRPIX2"]-=ymin
        h.pop("CHECKSUM",None);h.pop("DATASUM",None)
        hdus.append(fits.ImageHDU(data=arr,header=h))
    out=dest/filename
    with fits.HDUList(hdus) as fs:
        temp=out.with_suffix(".tmp");fs.writeto(temp,overwrite=True,checksum=True);temp.replace(out)
    return (*base,xmin,ymin,out.name,out.stat().st_size,hashlib.sha256(out.read_bytes()).hexdigest()),receipts


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("products",type=Path);p.add_argument("out",type=Path)
    p.add_argument("--planet",required=True);p.add_argument("--proposal",required=True)
    p.add_argument("--filters",nargs="+",required=True);p.add_argument("--half",type=int,default=256)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    with a.products.open() as f:
        rows=[r for r in csv.DictReader(f) if r['planet']==a.planet and r['proposal_id']==a.proposal
              and r['filters'] in a.filters and r['productFilename'].endswith('_cal.fits')]
    if not rows:raise ValueError("no products")
    with ThreadPoolExecutor(12) as pool:
        results=list(pool.map(lambda r:make_cut(r,a.half,a.out),rows))
    for name,fields,data in [
        ("cutouts","source filter visit midpoint_mjd_utc planet_x planet_y cutout_x0 cutout_y0 filename output_bytes output_sha256".split(),
         [r[0] for r in results]),
        ("ranges","source offset length sha256".split(),[receipt for r in results for receipt in r[1]])]:
        with (a.out/f"{name}.csv").open("w",newline="") as f:
            w=csv.writer(f,lineterminator="\n");w.writerow(fields);w.writerows(data)


if __name__=="__main__":
    main()
