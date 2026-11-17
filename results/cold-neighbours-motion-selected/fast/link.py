"""Blindly link bright W2-only epoch detections at 0.8–12 arcsec/yr; no CatWISE seed."""
import csv
import gzip
import hashlib
import io
import json
import math
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent
CACHE=Path(os.environ.get("UNTIMELY_CACHE",ROOT/"cache"))
BASE="https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/"
RADIUS=.4
REF=57170.
def get(name):
    p=CACHE/Path(name).name
    if not p.exists():
        with urlopen(BASE+name,timeout=150) as f:data=f.read()
        if not data.startswith(b"\x1f\x8b"):raise ValueError(name)
        p.write_bytes(data)
    data=p.read_bytes()
    return fits.getdata(io.BytesIO(gzip.decompress(data)),1),dict(url=BASE+name,sha256=hashlib.sha256(data).hexdigest())
def load(job):
    tile,a,d,obs=job;factor=math.cos(math.radians(d))
    epochs=[];refs=[]
    for o in obs:
        rows,rec=get(o["CATALOG_FILENAME"].removeprefix("./"));refs.append(rec)
        xy=np.column_stack([((rows["ra"]-a+180)%360-180)*factor*3600,(rows["dec"]-d)*3600])
        cut=(np.linalg.norm(xy,axis=1)<RADIUS*3600+150)&(rows["qf"]>.85)&(rows["dflux"]>0)&(rows["flux"]>=6*rows["dflux"])&(rows["flags_unwise"]==0)
        xy=xy[cut]; flux=rows["flux"][cut]; t=(o["MJDMEAN"]-REF)/365.25
        if len(xy):epochs.append((t,xy,cKDTree(xy),flux,rows["unwise_detid"][cut].astype(str)))
    epochs.sort(key=lambda e:e[0]);n=len(epochs);tracks=[]
    # Two seed pairs permit one early blended visit. No requirement that W1 be seen.
    start=[i for i,e in enumerate(epochs) if 56580<=e[0]*365.25+REF<=57400]
    end=[i for i,e in enumerate(epochs) if e[0]*365.25+REF>=58400]
    if len(start)<2 or len(end)<2:return [],dict(tile=tile,n_epochs=n,seeds=0,urls=refs)
    seeds=0
    for ia,ib in [(start[0],end[-1]),(start[1],end[-2])]:
        ta,x,tx,fx,_=epochs[ia];tb,y,ty,fy,_=epochs[ib];dt=tb-ta
        # Brightness cut uses native Vega nanomaggies (W2<=16 at both seed visits).
        brightflux=10**((22.5-16)/2.5)
        ix=np.where((fx>=brightflux)&(ty.query(x)[0]>2.0))[0]
        iy=np.where((fy>=brightflux)&(tx.query(y)[0]>2.0))[0]
        if not len(iy) or not len(ix):continue
        pairs=cKDTree(y[iy]).query_ball_point(x[ix],12*dt)
        for i,jlist in zip(ix,pairs):
            for k in jlist:
                j=iy[k];v=(y[j]-x[i])/dt;pos=x[i]-ta*v
                if np.linalg.norm(v)<.8 or np.linalg.norm(pos)>RADIUS*3600 or not .5<fx[i]/fy[j]<2:continue
                seeds+=1;sel=[]
                for seq,(t,xy,tr,flux,ids) in enumerate(epochs):
                    dist,num=tr.query(pos+t*v,distance_upper_bound=1.6)
                    if dist<1.6 and .4<flux[num]/fx[i]<2.5:sel.append((seq,num))
                if len(sel)<6:continue
                ts=np.array([epochs[b][0] for b,c in sel]);points=np.array([epochs[b][1][c] for b,c in sel])
                if np.ptp(ts)<5:continue
                co=np.linalg.lstsq(np.column_stack((np.ones(len(ts)),ts)),points,rcond=None)[0]
                resid=points-(co[0]+ts[:,None]*co[1]);rms=np.sqrt(np.mean(np.sum(resid**2,axis=1)))
                if rms>.9 or not .8<=np.linalg.norm(co[1])<=12 or np.linalg.norm(co[0])>RADIUS*3600:continue
                stationary=sum(e[2].query(co[0],distance_upper_bound=1.6)[0]<1.6 for e in epochs)
                if len(sel)-stationary<3:continue
                w2=22.5-2.5*math.log10(np.median([epochs[b][3][c] for b,c in sel]))
                tracks.append(dict(tile=tile,ra2015_deg=(a+co[0,0]/(3600*factor))%360,dec2015_deg=d+co[0,1]/3600,
                                   pmra=co[1,0],pmdec=co[1,1],w2=w2,n=len(sel),baseline=np.ptp(ts),rms=rms,
                                   fixed_epochs=stationary,ids=" ".join(epochs[b][4][c] for b,c in sel)))
    kept=[]
    for row in sorted(tracks,key=lambda r:(-r["n"],r["rms"])):
        if any(len(set(row["ids"].split())&set(s["ids"].split()))>=3 for s in kept):continue
        row["source_name"]="track"+tile+"_"+str(len(kept)+1);kept.append(row)
    print(tile,n,seeds,len(kept),flush=True)
    return kept,dict(tile=tile,ra=a,dec=d,radius_deg=RADIUS,n_epochs=n,seeds=seeds,urls=refs)
def main():
    CACHE.mkdir(exist_ok=True,parents=True)
    idx,rec=get("untimely_index-neo7.fits.gz")
    targets=sorted({r["source_id"][:8] for r in csv.DictReader((ROOT.parent/"pilot/shortlist.csv").open())})
    jobs=[]
    for tile in ["1342m076",*targets]:
        rows=idx[(idx["COADD_ID"].astype(str)==tile)&(idx["BAND"]==2)]
        bins=np.floor((rows["MJDMEAN"]-55197)/182.625)
        obs=[g[np.argmax(g["N_DET"])] for b in np.unique(bins) for g in [rows[bins==b]]]
        a,d=(133.7833,-7.2448) if tile=="1342m076" else (float(rows[0]["RA"]),float(rows[0]["DEC"]))
        jobs.append((tile,a,d,obs))
    allrows=[];meta=[]
    with ThreadPoolExecutor(3) as pool:
        for rr,m in pool.map(load,jobs):
            allrows.extend(rr);meta.append(m)
            cols=["source_name","tile","ra2015_deg","dec2015_deg","pmra","pmdec","w2","n","baseline","rms","fixed_epochs","ids"]
            with (ROOT/"tracks.csv").open("w") as f:
                w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(allrows)
            (ROOT/"manifest.json").write_text(json.dumps(dict(index=rec,fields=meta),indent=1)+"\n")
if __name__=="__main__":main()
