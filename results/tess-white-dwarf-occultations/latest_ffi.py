"""Read the latest 432 listed 200-s FFIs for the brightest overlapping clean target."""
import argparse
import csv
import hashlib
import io
import re
import json
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib.request import Request,urlopen
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from wd_search import search

BASE="https://stpubdata.s3.amazonaws.com/tess/public/ffi"
ROSTER="https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_110_ffic.sh"

def byte_range(url,start,end):
    for attempt in range(3):
        try:
            with urlopen(Request(url,headers={"Range":f"bytes={start}-{end}"}),timeout=25) as r:
                if r.status!=206 or not r.headers["Content-Range"].startswith(f"bytes {start}-{end}/"):
                    raise ValueError("range response mismatch")
                payload=r.read()
            if len(payload)!=end+1-start: raise ValueError("incomplete range")
            return payload
        except Exception:
            if attempt==2: raise
            time.sleep(attempt+1)

def one(name,star,keep_pixels=False):
    cam,ccd=star["camera"],star["ccd"]
    url=f"{BASE}/s0110/{name[4:8]}/{name[8:11]}/{cam}-{ccd}/{name}"
    rec=dict(filename=name)
    try:
        header=byte_range(url,0,25919)
        buf=io.BytesIO(header)
        ph=fits.Header.fromfile(buf)
        im=fits.Header.fromfile(buf)
        offset=buf.tell()
        assert im["BITPIX"]==-32 and im["CAMERA"]==int(cam) and im["CCD"]==int(ccd)
        assert im.get("SECTOR",int(name.split("-")[1][1:]))==110 and im.get("BSCALE",1)==1 and im.get("BZERO",0)==0
        w=WCS(im,naxis=2)
        x,y=map(float,w.all_world2pix(float(star["ra_obs"]),float(star["dec_obs"]),0))
        cx,cy=round(x),round(y)
        nx,ny=im["NAXIS1"],im["NAXIS2"]
        if not (4<=cx<nx-4 and 4<=cy<ny-4):
            if keep_pixels:
                return dict(rec, status="target off stampable image", x=x, y=y,
                            header_sha256=hashlib.sha256(header).hexdigest(),
                            quality=im["DQUALITY"], utc=im["DATE-OBS"],
                            time=(im["TSTART"]+im["TSTOP"])/2, flux=None)
            raise ValueError("edge of image")
        start=offset+4*((cy-4)*nx+cx-4)
        end=start+4*(8*nx+9)-1
        strip=byte_range(url,start,end)
        raw=np.frombuffer(strip,dtype=">f4")
        cube=raw[np.arange(9)[:,None]*nx+np.arange(9)]
        yy,xx=np.indices((9,9))
        rr=(xx+cx-4-x)**2+(yy+cy-4-y)**2
        ap=rr<=1.5**2;edge=(xx==0)|(xx==8)|(yy==0)|(yy==8)
        sky=float(np.nanmedian(cube[edge]))
        flux=float((cube[ap]-sky).sum())
        rec.update(start=start,end=end,header_sha256=hashlib.sha256(header).hexdigest(),strip_sha256=hashlib.sha256(strip).hexdigest(),
                   time=(im["TSTART"]+im["TSTOP"])/2,utc=im["DATE-OBS"],quality=im["DQUALITY"],x=x,y=y,
                   sky=sky,flux=flux if np.isfinite(flux) else None,aperture_pixels=int(ap.sum()))
        if keep_pixels:
            rec.update(stamp_origin=[cx-4,cy-4], stamp=cube.tolist(),
                       wcs=w.to_header(relax=True).tostring(sep="\n", endcard=False, padding=False))
    except Exception as exc:rec.update(error=str(exc))
    return rec

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("roster",type=Path);p.add_argument("footprint",type=Path);p.add_argument("measurements",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    if not a.roster.exists():a.roster.write_bytes(urlopen(ROSTER,timeout=100).read())
    tics={r["source_id"]:r["tic"] for r in json.loads(a.measurements.read_text()) if r["status"]=="included" and r["sector"]==109}
    star=min((r for r in csv.DictReader(a.footprint.open()) if r["GaiaEDR3"] in tics),key=lambda r:float(r["Gmag"]))
    pair=star["camera"]+"-"+star["ccd"]
    files=sorted(set(re.findall(r"-o (tess\d+-s0110-"+pair+r"-\d+-s_ffic\.fits)",a.roster.read_text())))
    if not files:raise SystemExit("missing detector roster")
    names=files[-432:]
    cache=a.out/"range_cache";cache.mkdir(exist_ok=True)
    def cached(n):
        path=cache/(n+".json")
        if path.exists():
            r=json.loads(path.read_text())
            if "error" not in r:return r
        r=one(n,star);path.write_text(json.dumps(r)+"\n");return r
    with ThreadPoolExecutor(9) as pool: rows=list(pool.map(cached,names))
    good=sorted([r for r in rows if not r.get("error") and r["quality"]==0 and r["flux"] is not None],key=lambda r:r["time"])
    (a.out/"latest_ffi_ranges.jsonl").write_text("".join(json.dumps(r)+"\n" for r in rows))
    summary=dict(tic=tics[star["GaiaEDR3"]],gaia=star["GaiaEDR3"],wd=star["WDJname"],roster_url=ROSTER,
                 roster_sha256=hashlib.sha256(a.roster.read_bytes()).hexdigest(),roster_ccd_files=len(files),
                 attempted=len(rows),read_ok=sum("error" not in r for r in rows),quality0=len(good),first_file=names[0],last_file=names[-1])
    if good:
        summary.update(first_good_utc=good[0]["utc"],last_good_utc=good[-1]["utc"],median_aperture_e_s=float(np.median([r["flux"] for r in good])),
                       search=search([r["time"] for r in good],[r["flux"] for r in good],200))
    (a.out/"latest_ffi_summary.json").write_text(json.dumps(summary,indent=1)+"\n");print(summary)
