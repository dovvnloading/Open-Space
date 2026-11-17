"""Independently extract epochal W2 detections near CatWISE motion tracks.

The full catalogs can be fetched again. The resulting table preserves every
unTimely detection within 25 arcsec of a selected source's 2015.405 position,
plus four 4-arcsec windows around the moving and rotated control predictions.
"""
import csv
import gzip
import hashlib
import io
import json
import math
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from time import sleep
from urllib.request import urlopen

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent
BASE = "https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/"
OUT = ROOT / "untimely"
CACHE = Path(os.environ.get("UNTIMELY_CACHE",OUT/"cache"))
IDX = "untimely_index-neo7.fits.gz"
KEEP = ["unwise_detid","ra","dec","MJDMEAN","EPOCH","FORWARD","flux","dflux",
        "dx","dy","qf","fracflux","rchi2","nm","primary","flags_unwise","flags_info"]

def fetch(name):
    cache = CACHE/Path(name).name
    if cache.exists(): return cache.read_bytes()
    for trial in range(3):
        try:
            with urlopen(BASE+name, timeout=120) as f: data=f.read()
            if not data.startswith(b"\x1f\x8b"): raise ValueError("expected gzip FITS")
            cache.write_bytes(data)
            return data
        except Exception:
            if trial==2: raise
            sleep(trial+1)

def digest(data):
    return hashlib.sha256(data).hexdigest()

def extract(job):
    name, sources = job
    data=fetch(name)
    with fits.open(io.BytesIO(gzip.decompress(data)), memmap=False) as hdus:
        table=hdus[1].data
        rows=[]
        for s in sources:
            ra0=float(s["ra2015_deg"]); dec0=float(s["dec2015_deg"])
            offa=(table["ra"]-ra0)*math.cos(math.radians(dec0))*3600
            offd=(table["dec"]-dec0)*3600
            mask=np.hypot(offa,offd) <= 25
            if len(table):
                dt=(float(table["MJDMEAN"][0])-57170)/365.25
                ma=dt*float(s["pmra"]); md=dt*float(s["pmdec"])
                if math.hypot(ma, md)>21:
                    for a,b in [(ma,md),(-md,ma),(-ma,-md),(md,-ma)]:
                        mask |= np.hypot(offa-a,offd-b)<=4
            for idx in np.flatnonzero(mask):
                r=table[idx]
                rows.append(dict(source=s["source_name"],**{k:r[k] for k in KEEP}))
    return dict(url=BASE+name,sha256=digest(data),length=len(data)),rows

def main():
    CACHE.mkdir(parents=True,exist_ok=True);OUT.mkdir(exist_ok=True)
    sources=list(csv.DictReader((ROOT/"shortlist.csv").open()))
    target_tiles={s["source_id"][:8] for s in sources if int(s["bright_tier"])}
    # The default replays all tiles with at least one primary color/motion
    # selection. Set UNTIMELY_EXTRA_TILES="" to replay just the bright tier.
    extra=os.environ.get("UNTIMELY_EXTRA_TILES","all")
    if extra=="all":
        target_tiles.update(s["source_id"][:8] for s in sources)
    else:
        target_tiles.update(filter(None,extra.split(",")))
    index=fetch(IDX); records=[dict(url=BASE+IDX,sha256=digest(index),length=len(index))]
    d=fits.getdata(io.BytesIO(gzip.decompress(index)),1)
    tiles=np.asarray(d["COADD_ID"]).astype(str)
    jobs=[]
    for tile in sorted(target_tiles):
        subs=d[(tiles==tile)&(d["BAND"]==2)]
        # Highest detection count for each disjoint 182.625-day slice.
        # Several epochs cover only a narrow strip near the ecliptic pole.
        bins=np.floor((subs["MJDMEAN"]-55197.)/182.625)
        src=[s for s in sources if s["source_id"][:8]==tile]
        for b in np.unique(bins):
            group=subs[bins==b]
            best=group[np.argmax(group["N_DET"])]
            jobs.append((best["CATALOG_FILENAME"].removeprefix("./"),src))
    found=[]
    with ThreadPoolExecutor(5) as ex:
        for met,rows in ex.map(extract,jobs):records.append(met);found.extend(rows)
    with (OUT/"detections.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=["source",*KEEP]);w.writeheader();w.writerows(found)
    (OUT/"manifest.json").write_text(json.dumps(records,indent=2)+"\n")
    print(f"{len(records)-1} epoch catalogs, {len(found)} detections, {len(target_tiles)} tiles")
if __name__=="__main__":main()
