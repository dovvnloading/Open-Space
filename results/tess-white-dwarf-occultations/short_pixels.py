"""Extract each released SPOC target pixel file for the Sector 109 cohort and search it."""
import argparse
import csv
import hashlib
import json
import shutil
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from wd_search import search
from ffi_acquire import propagated,tsv_rows
from ffi_photometry import deposit,gaia_t

def acquire(path,url):
    if not path.exists():
        tmp=path.with_suffix(".part")
        with urlopen(url,timeout=150) as r,tmp.open("wb") as f:
            shutil.copyfileobj(r,f,1024*1024)
        tmp.replace(path)
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def run(row,star,root,out):
    path=root/row["filename"]
    result=out/(path.stem+".json")
    if result.exists() and json.loads(result.read_text()).get("status") in {"searched","excluded_crowding"}: return
    url="https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/"+path.name
    rec=dict(**row,url=url)
    try:
        rec["sha256"]=acquire(path,url)
        rec["bytes"]=path.stat().st_size
        with fits.open(path) as h:
            assert int(h[0].header["TICID"])==int(row["tic"]) and int(h[0].header["SECTOR"])==int(row["sector"])
            d=h[1].data
            dt=float(h[1].header["TIMEDEL"])*86400
            assert abs(dt-int(row["cadence_s"]))<0.1
            wcs=WCS(h[2].header,naxis=2)
            t=np.asarray(d["TIME"],float)
            year=2000+(np.nanmedian(t)+2457000.-2451545.)/365.25
            sx,sy=map(float,wcs.world_to_pixel_values(*propagated(star,year)))
            yy,xx=np.indices(h[2].data.shape)
            ap=(xx-sx)**2+(yy-sy)**2<=1.5**2
            if ap.sum()<4 or sx<1.5 or sy<1.5 or sx>xx.max()-1.5 or sy>yy.max()-1.5: raise ValueError("truncated aperture")
            ns=tsv_rows(root/(star["GaiaEDR3"]+"-gaia.tsv"))
            if len([r for r in ns if r["Source"]==star["GaiaEDR3"]])!=1: raise ValueError("target Gaia row")
            fractions=[]
            for sigma in (.65,.85,1.05):
                base=deposit(xx[ap],yy[ap],sx,sy,sigma)
                for worst in (False,True):
                    others=0.
                    for n in ns:
                        if n["Source"]==star["GaiaEDR3"]:continue
                        nx,ny=wcs.world_to_pixel_values(*propagated(n,year))
                        others+=10**(-.4*(gaia_t(n,worst)-float(star["tmag"])-(.3 if worst else 0)))*deposit(xx[ap],yy[ap],nx,ny,sigma)/base
                    fractions.append(1/(1+others))
            rec.update(target_xy=[sx,sy],nominal_fraction=fractions[2],min_fraction=min(fractions),cadences=len(d),aperture_pixels=int(ap.sum()))
            if min(fractions)<.8:
                rec["status"]="excluded_crowding"
            else:
                # Calibrated SPOC target-pixel FLUX already has the pipeline sky subtracted.
                flux=np.asarray(d["FLUX"][:,ap],dtype=float).sum(axis=1)
                good=(d["QUALITY"]==0)&np.isfinite(t+flux)
                rec.update(start_btjd=float(min(t[good])),end_btjd=float(max(t[good])),status="searched",
                           median_e_s=float(np.median(flux[good])),search=search(t[good],flux[good],dt))
        print(path.name,rec["status"],flush=True)
    except Exception as exc:
        rec.update(status="error",error=repr(exc));print(path.name,rec["error"],flush=True)
    result.write_text(json.dumps(rec,indent=1)+"\n")

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("inventory",type=Path);p.add_argument("measurements",type=Path);p.add_argument("footprint",type=Path)
    p.add_argument("cache",type=Path);p.add_argument("out",type=Path);p.add_argument("--workers",type=int,default=3)
    a=p.parse_args();a.cache.mkdir(exist_ok=True,parents=True);a.out.mkdir(exist_ok=True,parents=True)
    sources={r["GaiaEDR3"]:r for r in csv.DictReader(a.footprint.open())}
    eligible={str(v["tic"]):dict(sources[v["source_id"]],tmag=v["tmag"]) for v in json.loads(a.measurements.read_text()) if v["status"]=="included" and v["sector"]==109}
    jobs=[r for r in csv.DictReader(a.inventory.open()) if r["tic"] in eligible]
    with ThreadPoolExecutor(a.workers) as p:
        list(p.map(lambda r:run(r,eligible[r["tic"]],a.cache,a.out),jobs))
