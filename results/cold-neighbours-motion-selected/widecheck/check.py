"""Compare wide-cap CatWISE tracks to independent unTimely W2 detections.

Use six fixed, disjoint time windows. Re-use optional cached native catalogs,
and otherwise stream and discard each complete catalog after local extraction.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import statistics
import tempfile
import shutil
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from time import sleep
from urllib.request import urlopen

import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parent
BASE="https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/"
INDEX="untimely_index-neo7.fits.gz"
WINDOWS=[(55197,55750),(56580,57024),(57024,57754),(57754,58466),(58466,58850),(58850,59216)]
FIELDS=["unwise_detid","ra","dec","MJDMEAN","EPOCH","FORWARD","flux","dflux","qf","fracflux","nm","flags_unwise"]
MODELS=["moving","stationary","rot90","rot180","rot270"]

def digest(d):return hashlib.sha256(d).hexdigest()
def fetch(path,cache):
    p=cache/path.rsplit("/",1)[-1] if cache is not None else None
    if p is not None and p.exists():return p.read_bytes()
    for trial in range(4):
        try:
            with urlopen(BASE+path,timeout=150) as r:raw=r.read()
            if not raw.startswith(b"\x1f\x8b"):raise ValueError(f"Invalid gzip resource {path}")
            return raw
        except Exception:
            if trial==3:raise
            sleep(2**trial)
def choose(t):
    out=[]
    for i,(lo,hi) in enumerate(WINDOWS):
        g=t[(t["MJDMEAN"]>=lo)&(t["MJDMEAN"]<hi)]
        if not len(g):continue
        # Prefer balanced scan directions when the requested pass covers >95%
        # of this tile. Rank the eligible pass by catalogue detection count.
        good=(g["FORWARD"]==i%2)&(g["FRAC_PIX_COV0"]<=.05)
        if not np.any(good):good=g["FRAC_PIX_COV0"]<=.05
        if np.any(good):g=g[good]
        out.append(g[np.argmax(g["N_DET"])])
    return out
def tile_job(tile,src,epochs,cache,state):
    ids=sorted(s["source_id"] for s in src)
    saved=state/(tile+".json")
    signature=digest(json.dumps(src,sort_keys=True).encode())
    if saved.exists():
        old=json.loads(saved.read_text())
        if old.get("signature")==signature:return old
    meta=[];detections=[]
    for ep in epochs:
        path=ep["CATALOG_FILENAME"].removeprefix("./")
        data=fetch(path,cache)
        meta.append(dict(url=BASE+path,sha256=digest(data),length=len(data),mjd=round(float(ep["MJDMEAN"]),5),forward=int(ep["FORWARD"]),tile_fraction_zero_cov=round(float(ep["FRAC_PIX_COV0"]),5)))
        # Crowded Galactic-plane tiles can exceed 200 MB decompressed. Use
        # a temporary memory-mapped FITS file instead of full byte copies in
        # each worker. The mapping stays usable after the temporary unlink.
        with tempfile.NamedTemporaryFile(dir=state) as uncompressed:
            with gzip.GzipFile(fileobj=io.BytesIO(data)) as stream:
                shutil.copyfileobj(stream,uncompressed,length=1024*1024)
            uncompressed.flush()
            table=fits.getdata(uncompressed.name,1,memmap=True)
        for s in src:
            dec=float(s["dec2015_deg"]);ra=float(s["ra2015_deg"])
            xx=((table["ra"]-ra+180)%360-180)*3600*math.cos(math.radians(dec))
            yy=(table["dec"]-dec)*3600
            dt=(float(ep["MJDMEAN"])-57170)/365.25
            a,b=dt*float(s["pmra"]),dt*float(s["pmdec"])
            mask=np.hypot(xx,yy)<=25
            if math.hypot(a,b)>21:
                for x,y in [(a,b),(-b,a),(-a,-b),(b,-a)]:mask|=np.hypot(xx-x,yy-y)<=4
            # Keep local source rows relevant for either assessment of depth
            # or association with a trial track.
            mask&=(table["qf"]>.85)&(table["dflux"]>0)
            for k in np.flatnonzero(mask):
                r=table[k]
                row={col:(str(r[col]) if col=="unwise_detid" else r[col].item()) for col in FIELDS}
                detections.append(dict(source_name=s["source_name"],**row))
        del table,data
    result=dict(tile=tile,signature=signature,source_ids=ids,epochs=meta,detections=detections)
    tmp=saved.with_suffix(".tmp");tmp.write_text(json.dumps(result,separators=(",",":")));tmp.replace(saved)
    return result
def score(s,rows,n):
    ra=float(s["ra2015_deg"]);dec=float(s["dec2015_deg"]);ma=float(s["pmra"]);md=float(s["pmdec"])
    velocities=[(ma,md),(0,0),(-md,ma),(-ma,-md),(md,-ma)]
    matches={m:{} for m in MODELS};noise=defaultdict(list)
    for d in rows:
        x=((d["ra"]-ra+180)%360-180)*3600*math.cos(math.radians(dec));y=(d["dec"]-dec)*3600
        mjd=d["MJDMEAN"];t=(mjd-57170)/365.25
        if d["nm"]>=3 and math.hypot(x,y)<=25:noise[mjd].append(22.5-2.5*math.log10(4*d["dflux"]))
        if d["flux"]/d["dflux"]<4:continue
        for model,(a,b) in zip(MODELS,velocities):
            dist=math.hypot(x-t*a,y-t*b)
            if dist<=1.5 and dist<matches[model].get(mjd,(99,None))[0]:
                matches[model][mjd]=(dist,d)
    ds=[statistics.median(v) for v in noise.values()];mv=matches["moving"];m=len(mv)
    rival=max(len(matches[k]) for k in MODELS[1:])
    times=sorted(mv); late=sum(t>58465 for t in times)
    span=(times[-1]-times[0])/365.25 if times else 0
    out=dict(field=s["field"],source_name=s["source_name"],source_id=s["source_id"],w2=s["w2"],
        w1w2=s["w1w2"],mu=s["mu"],n_files=n,local_depth_epochs=len(ds),
        local_four_sigma_median=round(statistics.median(ds),3) if ds else "",
        depth_reaches_w2=sum(v>=float(s["w2"]) for v in ds),
        baseline_yr=round(span,3),post_catwise=late,distinct_directions=len({t[1]["FORWARD"] for t in mv.values()}),
        moving_advantage=m-rival,priority=int(m>=3 and span>=5 and late>=1 and m-rival>=2),
        moving_ids=" ".join(mv[t][1]["unwise_detid"] for t in times))
    out.update({k:len(v) for k,v in matches.items()})
    # The nominal CatWISE reference position need not coincide with a
    # contaminant. Compare free reference-point shifts on the same grid for
    # the reported velocity and a stationary source (descriptive scores).
    if m>=3:
        grid=np.arange(-3,3.01,.375)
        xg,yg=np.meshgrid(grid,grid);inside=xg*xg+yg*yg<=9
        gx=xg[inside];gy=yg[inside]
        for tag,(va,vd) in [("moving",(ma,md)),("stationary",(0,0))]:
            by_epoch={}
            for d in rows:
                if d["flux"]/d["dflux"]<4:continue
                mjd=d["MJDMEAN"];dt=(mjd-57170)/365.25
                ox=((d["ra"]-ra+180)%360-180)*3600*math.cos(math.radians(dec))-va*dt
                oy=(d["dec"]-dec)*3600-vd*dt
                if abs(ox)>4.5 or abs(oy)>4.5:continue
                hits=(gx-ox)**2+(gy-oy)**2<=2.25
                if mjd in by_epoch:by_epoch[mjd]|=hits
                else:by_epoch[mjd]=hits
            out["offset_"+tag]=int(np.max(np.sum(list(by_epoch.values()),axis=0))) if by_epoch else 0
    else:
        out.update(offset_moving="",offset_stationary="")
    # Apply the same thresholds after exchanging the trial motion and each
    # rotated control; the stationary path remains a shared competitor.
    for ctrl in ("rot90","rot180","rot270"):
        hits=sorted(matches[ctrl]);best=max(len(v) for k,v in matches.items() if k!=ctrl)
        out[ctrl+"_passes"]=int(len(hits)>=3 and (hits[-1]-hits[0])/365.25>=5
            and any(t>58465 for t in hits) and len(hits)-best>=2)
    return out
def write_csv(path,rows):
    if not rows:return
    temp=path.with_suffix(".tmp")
    with temp.open("w") as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    temp.replace(path)
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--catalog",type=Path,default=ROOT.parent/"wide"/"selected.csv")
    p.add_argument("--literature",type=Path,default=ROOT.parent/"catalogues"/"wide_matches.csv")
    p.add_argument("--cache",type=Path,default=os.environ.get("UNTIMELY_CACHE"))
    p.add_argument("--state",type=Path,default=ROOT/"per_tile")
    p.add_argument("--workers",type=int,default=3)
    a=p.parse_args();a.state.mkdir(exist_ok=True,parents=True)
    input_bytes=a.catalog.read_bytes();sources=list(csv.DictReader(io.StringIO(input_bytes.decode())))
    known={s["source_name"] for s in csv.DictReader(a.literature.open())}
    index=fetch(INDEX,a.cache);hash_index=digest(index);tbl=fits.getdata(io.BytesIO(gzip.decompress(index)),1)
    tiles=np.asarray(tbl["COADD_ID"]).astype(str)
    groups=defaultdict(list)
    for s in sources:groups[s["source_id"][:8]].append(s)
    extracted=[]
    with ThreadPoolExecutor(a.workers) as ex:
        todo={ex.submit(tile_job,k,v,choose(tbl[(tiles==k)&(tbl["BAND"]==2)]),a.cache,a.state):k for k,v in groups.items()}
        for fut in as_completed(todo):extracted.append(fut.result())
    records=[];allrows=[];summ=[]
    for tile in sorted(extracted,key=lambda x:x["tile"]):
        records.extend(tile["epochs"]);allrows.extend(tile["detections"])
        d=defaultdict(list)
        for x in tile["detections"]:d[x["source_name"]].append(x)
        for s in groups[tile["tile"]]:
            v=score(s,d[s["source_name"]],len(tile["epochs"]))
            v["listed_match"]=int(s["source_name"] in known);summ.append(v)
    write_csv(ROOT/"summary.csv",sorted(summ,key=lambda s:s["source_name"]))
    # Preserve the exact CSV bytes, then compress with a deterministic header.
    payload=io.StringIO(newline="")
    writer=csv.DictWriter(payload,fieldnames=allrows[0]);writer.writeheader();writer.writerows(allrows)
    (ROOT/"detections.csv.gz").write_bytes(gzip.compress(payload.getvalue().encode(),compresslevel=6,mtime=0))
    info=dict(catalog_sha256=digest(input_bytes),literature_sha256=digest(a.literature.read_bytes()),index_url=BASE+INDEX,index_sha256=hash_index,catalogues=records)
    (ROOT/"manifest.json").write_text(json.dumps(info,indent=1)+"\n")
    print(len(sources),len(groups),len(records),len(allrows),sum(v["priority"] for v in summ),flush=True)
if __name__=="__main__":main()
