"""Screen released SPOC simple-aperture light curves for deep short dips.

The full 20/120-s product roster is scanned. The separate Sector 109 cohort also
has an independent target-pixel extraction. This extension has no rate estimator.
"""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen
from time import sleep
import numpy as np
from astropy.io import fits
from wd_search import mad, flatten, box_scores, event_seeds, period_seeds, periodic_candidates

def scan(time,flux,cadence):
    dt=cadence/86400
    trend=flatten(time,flux)
    good=np.isfinite(trend)
    t=time[good];f=flux[good]/trend[good]-1
    if len(f)<100 or np.ptp(t) < .16: return [],[]
    noise=mad(f);periodic=[];isolated=[]
    if not np.isfinite(noise) or noise<=0:return [],[]
    for minutes in (2,4,8,16,32):
        width=max(minutes/1440,2*dt)
        depth,z=box_scores(t,f,dt,width,noise)
        good=(z>=4)&(depth>=.02)
        # This shallow-stage gate requires 3 individually significant eclipses.
        if not good.any():continue
        ix=event_seeds(t,np.where(good,z,np.nan),width)
        for i in ix:
            if z[i]>=10 and depth[i]>=.05:
                isolated.append(dict(epoch_btjd=float(t[i]),duration_min=width*1440,depth=float(depth[i]),score=float(z[i])))
        if len(ix)>=3:
            use=ix[:18]
            trial=period_seeds(t[use],z[use],width,np.ptp(t))
            for r in periodic_candidates(t,depth,z,dt,width,trial):
                if r["depth"]>=.02 and r["score"]>=10:periodic.append(r)
    return periodic,isolated

def product(row,exact,out):
    name=row["filename"].replace("_tp.fits","_lc.fits").replace("_fast-tp.fits","_fast-lc.fits")
    path=out/(name+".json")
    if path.exists() and "error" not in json.loads(path.read_text()):return
    tic=row["tic"].zfill(16)
    url=f"https://stpubdata.s3.amazonaws.com/tess/public/tid/s{int(row['sector']):04d}/"+"/".join(tic[i:i+4] for i in (0,4,8,12))+"/"+name
    r=dict(tic=row["tic"],sector=int(row["sector"]),cadence_s=int(row["cadence_s"]),filename=name,url=url)
    if row["tic"] not in exact:
        r["status"]="position-only TIC match; excluded"
    else:
        try:
            for attempt in range(3):
                try:
                    with urlopen(url,timeout=60) as resp: raw=resp.read()
                    break
                except OSError:
                    if attempt==2:raise
                    sleep(attempt+1)
            r.update(bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
            with fits.open(io.BytesIO(raw),memmap=False) as h:
                hdr=h[1].header
                if int(h[0].header["TICID"])!=int(tic) or int(h[0].header["SECTOR"])!=r["sector"] or abs(86400*hdr["TIMEDEL"]-r["cadence_s"])>1:raise ValueError("product identity mismatch")
                frac=hdr.get("CROWDSAP")
                r["crowdsap"]=frac
                if frac is None or not .8<=frac<=1.:r["status"]="excluded by SPOC crowding"
                else:
                    d=h[1].data;t=d["TIME"];f=d["SAP_FLUX"]
                    use=(d["QUALITY"]==0)&np.isfinite(t+f)
                    t,ix=np.unique(t[use],return_index=True);f=f[use][ix].astype(float)
                    if len(t)<100 or np.median(f)<=0:r["status"]="insufficient flux"
                    else:
                        p,i=scan(t,f,r["cadence_s"]);r.update(status="searched",quality0=len(t),periodic=p,isolated=i)
        except Exception as exc:r["error"]=str(exc)
    path.write_text(json.dumps(r,allow_nan=False)+"\n")

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("inventory",type=Path);p.add_argument("crossmatch",type=Path);p.add_argument("out",type=Path);p.add_argument("--workers",type=int,default=5);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    exact={r["tic"] for r in csv.DictReader(a.crossmatch.open()) if r["status"]=="source_id_and_position"}
    rows=list(csv.DictReader(a.inventory.open()))
    with ThreadPoolExecutor(a.workers) as ex:list(ex.map(lambda r:product(r,exact,a.out),rows))
