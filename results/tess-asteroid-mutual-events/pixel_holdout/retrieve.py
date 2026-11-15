"""Reconstruct the moving-target stamps using calibrated MAST full-frame images on AWS.

Use a manifest containing the exact file name of each accepted FFI. Request the
FITS header plus one contiguous byte range for thirteen rows around the target.
MPC C57 astrometry supplies the time-dependent sky coordinate on the detector.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import gzip
import io
from pathlib import Path
from urllib.request import Request, urlopen

from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import numpy as np

ROOT=Path(__file__).resolve().parent
BUCKET="https://stpubdata.s3.amazonaws.com/tess/public/ffi"

def byte_range(url,start,length):
    request=Request(url,headers={"Range":f"bytes={start}-{start+length-1}"})
    for attempt in range(3):
        try:
            with urlopen(request,timeout=45) as response:
                if response.status!=206:
                    raise ValueError(f"Range request returned {response.status}")
                return response.read()
        except OSError:
            if attempt==2:
                raise

def checksum(path):
    h=hashlib.sha256()
    with np.load(path) as values:
        for key in ("t","x","y","q","flux"):
            data=np.asarray(values[key],dtype="<i8" if key=="q" else "<f8")
            h.update(data.tobytes(order="C"))
    return h.hexdigest()

def reconstruct(key,ephemeris,destination):
    filename=key.removesuffix(".npz")
    name=filename.split("-")
    stamp=name[0][4:]
    url=f"{BUCKET}/{name[1]}/{stamp[:4]}/{stamp[4:7]}/{name[2]}-{name[3]}/{filename}"
    with io.BytesIO(byte_range(url,0,23040)) as stream:
        fits.Header.fromfile(stream)
        header=fits.Header.fromfile(stream)
        offset=stream.tell()
    t=(Time(header["DATE-OBS"],scale="utc").jd+Time(header["DATE-END"],scale="utc").jd)/2
    times,ra,dec=ephemeris.T
    ra_deg,dec_deg=np.interp(t,times,ra),np.interp(t,times,dec)
    col,row=WCS(header,naxis=2).all_world2pix(ra_deg,dec_deg,0)
    x,y=int(np.rint(col)),int(np.rint(row))
    width=header["NAXIS1"]
    if header["BITPIX"]!=-32 or not (6<=x<width-6 and 6<=y<header["NAXIS2"]-6):
        raise ValueError("Unexpected FFI geometry")
    start=offset+4*((y-6)*width+x-6)
    block=byte_range(url,start,4*(12*width+13))
    pix=np.stack([np.frombuffer(block,dtype=">f4",count=13,offset=r*width*4) for r in range(13)]).astype(float)
    np.savez_compressed(destination/(filename+".npz"),t=t,x=col-x,y=row-y,q=int(header.get("DQUALITY",0)),flux=pix)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir",type=Path)
    parser.add_argument("--workers",type=int,default=6)
    args=parser.parse_args()
    with gzip.open(ROOT/"ffi_manifest.csv.gz","rt") as f:
        frames=list(csv.DictReader(f))
    with gzip.open(ROOT/"mpc_ephemeris.csv.gz","rt") as f:
        astrometry=list(csv.DictReader(f))
    for visit in sorted({r["visit"] for r in frames}):
        asteroid=visit.split("_")[0]
        eph=np.array([[float(r[k]) for k in ("jd_utc","ra_deg","dec_deg")] for r in astrometry if r["mpc"]==asteroid])
        dest=args.workdir/("stamps"+visit)
        dest.mkdir(exist_ok=True,parents=True)
        def fetch_and_verify(row):
            path=dest/(row["ffi"]+".npz")
            if not path.exists():
                reconstruct(row["ffi"],eph,dest)
            if checksum(path)!=row["sha256_numeric_payload"]:
                raise ValueError(f"Numerical stamp checksum changed for {path.name}")
        selected=[r for r in frames if r["visit"]==visit]
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for _ in pool.map(fetch_and_verify,selected):
                pass
if __name__=="__main__":
    main()
