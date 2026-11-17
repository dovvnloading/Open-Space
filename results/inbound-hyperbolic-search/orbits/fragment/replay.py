"""Separate demonstrable pre-fragmentation positions of comet C/1996 J1."""
import base64
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import html
import json
from pathlib import Path
import re
import sys
import urllib.request

import fitz
from jplephem.spk import SPK
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/"slow"))
from replay import acquire, conic_seed, fit, KERNEL_SHA, PLANETS, STATIONS
CACHE=ROOT/"cache"
OUT=ROOT/"evidence"
CACHE.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)
acquire.CACHE=fit.CACHE=CACHE
fit.DATA=OUT
CUTOFF="1996 10 01"
KEYS={"J96J010","J96J01a","J96J01b"}
MPEC="https://www.minorplanetcenter.net/mpec/J96/J96K02.html"
MONTHS=["19960601","19960701","19960730","19960828","19960927"]
PDF="https://www.minorplanetcenter.net/iau/ECS/MPCArchive/1996/MPC_{}.pdf"

def digest(b): return hashlib.sha256(b).hexdigest()

def ident(line):
    return (line[15:32].strip(),line[32:44].strip(),line[44:56].strip(),line[77:80])

def squeeze(s):
    return "".join(s.replace("−","-").replace("–","-").split())

def download_mpc(url):
    record={"url":url}
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"Open-Space comet research/1"}),timeout=22) as r:
            b=r.read()
            record.update(final_url=r.url,http_date=r.headers.get("Date"),sha256=digest(b),bytes=len(b))
        if url.endswith(".pdf"):
            with fitz.open(stream=b,filetype="pdf") as doc:
                pages=[]
                for i,p in enumerate(doc):
                    words=p.get_text("words")
                    # Reconstruct printed rows by baseline. Two-column pages
                    # can print two observation sequences at different x's.
                    by_y={}
                    for x0,y0,x1,y1,word,*_ in words:
                        by_y.setdefault(round((y0+y1)/3),[]).append((x0,word))
                    for vals in by_y.values():
                        line=squeeze(" ".join(v[1] for v in sorted(vals)))
                        if "1996J1" in line:
                            pages.append((i+1,line))
                return record,pages
        clean=html.unescape(re.sub("<[^>]*>","",b.decode("latin1")))
        lines=[s[:80] for s in clean.splitlines() if len(s)>=80 and s[5:12]=="J96J010"]
        record["original_designation_rows"]=len(lines)
        return record,lines
    except Exception as exc:
        record["failed"]=str(exc)
        return record,[]

def model(key,lines,el,eph):
    raw=("\n".join(lines)+"\n").encode("ascii")
    (OUT/(key+".obs")).write_bytes(raw)
    ans={"sha256":digest(raw),"dates":[min(x[15:32].strip() for x in lines),max(x[15:32].strip() for x in lines)],
         "lines":len(lines),"models":{}}
    try:
        t,observer,east,north,sig,groups,dropped=fit.observations(key,eph)
        epoch=float(np.ceil(t.max())+1)
        eph_short=fit.ephemeris(eph,t.min()-fit.OBS_TIME_LIMIT-2,epoch+1,.25)
        seed=conic_seed(el,el["epoch"],epoch,eph)
        ans.update(epoch_jd_tdb=epoch,supported=len(t),dropped=dropped,arc_days=float(np.ptp(t)),
                   stations=len({g[0] for g in groups}),station_nights=len(set(groups)))
        for mode,ng in [("gravity",0),("radial_r2",1),("marsden_outward",2)]:
            r={}
            ans["models"][mode]=r
            try:
                x,cov,mask,res,red,pred=fit.fit(np.r_[seed,np.zeros(ng)],epoch,t,observer,east,north,sig,eph_short,
                                           outward=(ng==2),return_predictor=True)
                if ng==0: seed=x.copy()
                r.update(used=int(mask.sum()),rms_coordinate_arcsec=float(np.sqrt(np.mean(res[mask]**2))),
                         reduced_chi2=float(red),state_ssb_icrf=x.tolist(),covariance=cov.tolist())
                r["boundaries"]=fit.inbound(key,x,cov,epoch,eph)
                if ng==0:
                    r["holdouts"]=fit.check_holdouts(key,x,epoch,eph,t,observer,east,north,sig,groups,eph_short)
            except Exception as exc: r["failed"]=str(exc)
    except Exception as exc: ans["failed"]=str(exc)
    return ans

def main(src):
    receipt=json.loads((src/"source.json").read_text())
    if digest((src/"CmtObs.txt.gz").read_bytes())!=receipt["sha256_compressed"]:
        raise ValueError("snapshot checksum mismatch")
    records=[]
    with gzip.open(src/"CmtObs.txt.gz","rt",encoding="ascii",errors="replace") as f:
        for line in f:
            if line[5:12] in KEYS and len(line)>=80 and line[15:25]<CUTOFF and line[15:19]=="1996":
                records.append(line[:80])
    with ThreadPoolExecutor(max_workers=4) as pool:
        downloads=list(pool.map(download_mpc,[MPEC]+[PDF.format(v) for v in MONTHS]))
    paper=downloads[0][1]
    original={ident(l):l for l in paper}
    evidence={str(ident(l)):["MPEC 1996-K02"] for l in paper}
    all_pre={ident(l):l for l in records}
    all_pre.update(original)
    for line in records:
        obs=ident(line)
        if line[5:12]=="J96J010":
            original[obs]=line
            evidence.setdefault(str(obs),[]).append("bulk parent designation")
        stamp=squeeze("".join(obs[:3]))
        for rec,pages in downloads[1:]:
            for page,printed in pages:
                at=printed.find(stamp)
                if at<0 or "1996J1" not in printed[max(0,at-20):at] or obs[-1] not in printed[at+len(stamp):at+len(stamp)+18]:
                    continue
                evidence.setdefault(str(obs),[]).append(rec["url"]+"#page="+str(page))
                original[obs]=line
    queries=[]
    if digest(acquire.retrieve(PLANETS,"de440.bsp",queries))!=KERNEL_SHA: raise ValueError("DE440 hash mismatch")
    if b"Greenwich" not in acquire.retrieve(STATIONS,"stations.txt",queries): raise ValueError("bad observatory table")
    data=json.loads(gzip.decompress((ROOT.parents[1]/"catalog/sbdb.json.gz").read_bytes()))
    row=next(dict(zip(data["fields"],l)) for l in data["data"] if l[data["fields"].index("pdes")]=="1996 J1-B")
    el={k:float(row[k]) for k in ["q","e","i","om","w","tp","epoch"]}
    out={"cutoff_exclusive":CUTOFF,"split_reference":"https://ntrs.nasa.gov/citations/20060034802",
         "snapshot":receipt,"source_receipts":queries+[r for r,v in downloads],"pre_cutoff_packed_counts":dict(Counter(l[5:12] for l in records)),
         "crossmatches":evidence,"sets":{},"independent_conic_seed_spkid":row["spkid"]}
    with SPK.open(str(CACHE/"de440.bsp")) as eph:
        for name,data in [("documented_parent",original),("all_pre_cutoff_sensitivity",all_pre)]:
            selected=sorted(data.values(),key=lambda l:ident(l))
            if set(original)==set(all_pre) and name.startswith("all_pre"):
                out["sets"][name]={"identical_to":"documented_parent"}
                continue
            r=model(name,selected,el,eph) if selected else {"failed":"no matching records"}
            r["mpc80"]=selected
            out["sets"][name]=r
    raw=(json.dumps(out,indent=2,allow_nan=False)+"\n").encode()
    (OUT/"results.json").write_bytes(raw)
    blob=base64.b64encode(gzip.compress(raw,mtime=0)).decode()
    for i in range(0,len(blob),4500):
        print("PARENT_DATA",i,blob[i:i+4500],flush=True)

if __name__=="__main__": main(Path(sys.argv[1]))
