"""Uncapped second-pass spectroscopy of eligible, previously unresolved HETDEX hypotheses."""
import argparse
from collections import Counter, defaultdict, deque
from contextlib import ExitStack
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
from hashlib import sha256
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import tempfile
import time

import numpy as np
import requests
from astropy.io import fits
from scipy.spatial import cKDTree

from desi_lines import TAP, SPARCL, measure_arrays
from vetting import atom, chord, number, separation

RADIUS=1.25
BATCH=32
METALS=("civ1","civ2","oiii_uv1","oiii_uv2","ciii1","ciii2")
VIZ=("https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync",
     "https://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync")
NER="https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix"
READ=(5,19)


def dump(obj):
    return json.dumps(atom(obj),allow_nan=False,separators=(",",":"))


def lines(path):
    with gzip.open(path,"rt") as f:
        for line in f:
            yield json.loads(line)


def row_key(p):
    return (p["row_index"],p["channel"],p.get("stream") or "",float(p["wave"]).hex())


def row_z(p):
    return number(p.get("z_heii" if p["channel"]=="heii" else "z_lya"))


def strong_repeat(p,v):
    shot=p.get("shotid")
    prim={r["shotid"] for r in v.get("repeats",[]) if number(r.get("snr"))>=4 and r.get("shotid")!=shot}
    sec={r["shotid"] for r in v.get("secondary_repeats",[]) if number(r.get("snr"))>=4 and r.get("shotid")!=shot}
    return sorted(prim & sec)


def eligible(p,v):
    if v["disposition"]!="unresolved" or p["channel"] not in ("ew","heii"):
        return None
    if "shuffl" in str(p.get("stream","")).lower() or number(p.get("pair_offset_A"),0)!=0:
        return None
    f=number(p.get("flux"))
    if f<=0 or p.get("single_pixel") or p.get("sky_near") or p.get("badcal"):
        return None
    if p["channel"]=="ew":
        imaging=number(p.get("imaging_rank_A"),0)
        im=(v.get("imaging_EW") or {}).get("ew0_lower_A")
        if max(number(p.get("ew0_lower"),0),imaging,number(im,0))<=240:
            return None
        return dict(ew0_lower_A=p.get("ew0_lower"),imaging_ew0_lower_A=max(imaging,number(im,0)))
    partners=[v for v in (p.get("pair_lya"),p.get("lya_velocity_best")) if
              v and number(v.get("snr"))>=5 and number(v.get("flux"))>0]
    if not partners:
        return None
    # Prefer the directly fitted velocity/width if present.
    ly=partners[0]
    ratio=f/number(ly["flux"])
    width=2.355*number(p.get("sigma_A"),math.inf)*299792.458/p["wave"]
    if ratio<.1 or not 0<width<900 or not 1.9<=row_z(p)<=2.4:
        return None
    fam=p.get("heii_metal_components") or p.get("heii_family") or p.get("pair_metals") or {}
    if any(number(fam.get(k,{}).get("snr"))>=3 for k in METALS):
        return None
    coverage={}
    for ion,keys in (("CIV",METALS[:2]),("OIII]",METALS[2:4]),("CIII]",METALS[4:])):
        m=[fam.get(k,{}) for k in keys]
        ok=all(r.get("status")=="measured" and number(r.get("limit3"))>0 for r in m)
        ratio_limit=sum(r["limit3"] for r in m)/f if ok else None
        status="measured" if ok else "/".join(r.get("status","unreported") for r in m)
        coverage[ion]=dict(status=status,component_sum_upper3_over_heii=ratio_limit)
    return dict(heii_lya=ratio,lya_snr=ly["snr"],fwhm_kms=width,metals=coverage)


def select(indir):
    found=[]
    counters=Counter()
    for line,row in enumerate(lines(indir/"vetted.jsonl.gz"),1):
        p,v=row["peak"],row["vet"]
        counters["ledger_rows"]+=1
        criteria=eligible(p,v)
        if criteria is None:
            continue
        counters[p["channel"]+"_eligible"]+=1
        found.append(dict(first_pass_line=line,peak=p,first=v,criteria=criteria,
                     first_pass_external="SPARCL" in v,same_visit_two_line_repeats=strong_repeat(p,v)))
    if not found:
        return found,counters
    # Stable group id is for cross-release/visit review only. Each constituent
    # coordinate is queried and each line/row keeps its own adjudication key.
    xyz=chord([r["peak"]["ra"] for r in found],[r["peak"]["dec"] for r in found])
    tree=cKDTree(xyz)
    groups=np.arange(len(found))
    for i,j in tree.query_pairs(2*np.sin(np.deg2rad(1./3600)/2)):
        a,b=found[i]["peak"],found[j]["peak"]
        if a["channel"]!=b["channel"] or abs(a["wave"]-b["wave"])>3.:
            continue
        ii,jj=int(i),int(j)
        while ii!=groups[ii]:
            ii=groups[ii]
        while jj!=groups[jj]:
            jj=groups[jj]
        groups[max(ii,jj)]=min(ii,jj)
    for i,r in enumerate(found):
        j=i
        while j!=groups[j]:
            j=groups[j]
        r["group"]=int(j)
        r["key"]=row_key(r["peak"])
    return found,counters


def query_catalog(points,service):
    field="ra,dec" if service=="SPARCL" else 'RAICRS,DEICRS'
    where=" OR ".join(f"1=CONTAINS(POINT('ICRS',{field}),CIRCLE('ICRS',{ra:.9f},{dec:.9f},{RADIUS/3600:.11f}))"
                      for ra,dec in points)
    if service=="SPARCL":
        prefix="SELECT sparcl_id,ra,dec,data_release,specid,spectype,redshift,redshift_warning FROM sparcl.main"
        urls=[(TAP,"post"),(TAP,"get")]
        required={"sparcl_id","redshift","ra","dec","specid"}
    else:
        prefix=('SELECT TargetID, RAICRS, DEICRS, Survey, Prog, Healpix, z, ZWARN, OType, delChi2 '
                'FROM "V/161/zcatdr1"')
        urls=[(VIZ[0],"post"),(VIZ[1],"get")]
        required={"targetid","raicrs","deicrs","healpix","z","zwarn"}
    query=prefix+" WHERE ("+where+")"
    params={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","MAXREC":"10000","QUERY":query}
    receipt=dict(service=service,query=query,attempts=[])
    for url,method in urls:
        try:
            response=(requests.post(url,data=params,timeout=READ) if method=="post" else
                      requests.get(url,params=params,timeout=READ))
            response.raise_for_status()
            table=csv.DictReader(io.StringIO(response.content.decode("utf-8-sig")))
            cols={n.strip().lower() for n in table.fieldnames or []}
            if not required.issubset(cols):
                raise ValueError("Expected catalogue header unavailable: "+response.text[:165])
            rows=[{k.strip().lower():v for k,v in r.items()} for r in table]
            if len(rows)>=10000:
                raise ValueError("TAP row limit reached")
            if service=="CDS_DESI_DR1":
                for r in rows:
                    r.update(ra=r["raicrs"],dec=r["deicrs"],data_release="DESI-DR1",specid=r["targetid"],
                             spectype=r["otype"],redshift=r["z"],redshift_warning=r["zwarn"])
            receipt["attempts"].append(dict(url=url,method=method,status="complete",bytes=len(response.content),
                                         sha256=sha256(response.content).hexdigest(),rows=len(rows)))
            return rows,receipt
        except (OSError,ValueError,TypeError,csv.Error) as exc:
            receipt["attempts"].append(dict(url=url,method=method,status="failed",error=str(exc)[:280]))
            time.sleep(.15)
    return None,receipt


def coadd_url(m):
    survey=m["survey"].lower().strip()
    prog=m["prog"].lower().strip()
    pixel=int(m["healpix"])
    if not re.fullmatch("[a-z0-9]+",survey+prog) or not 0<=pixel<=49151:
        raise ValueError("Invalid DESI coadd identity")
    return f"{NER}/{survey}/{prog}/{pixel//100}/{pixel}/coadd-{survey}-{prog}-{pixel}.fits"


class RangeFile(io.RawIOBase):
    """Seekable FITS stream retaining downloaded bytes and range hashes."""
    def __init__(self,url):
        super().__init__()
        self.url=url
        self.session=requests.Session()
        self.pos=0
        self.deadline=time.monotonic()+110
        self.size=None
        self.parts={}
        self.receipts=[]
        self.fetch(0,2879)

    def fetch(self,start,end):
        for attempt in range(2):
            try:
                remain=self.deadline-time.monotonic()
                if remain<1:
                    raise TimeoutError("DESI range-read budget exhausted")
                r=self.session.get(self.url,headers={"Range":f"bytes={start}-{end}","Accept-Encoding":"identity"},
                                   timeout=(min(READ[0],remain),min(READ[1],remain)),stream=True)
                with r:
                    if r.status_code!=206:
                        raise ValueError(f"HTTP Range unavailable ({r.status_code})")
                    size=int(r.headers["Content-Range"].split("/")[-1])
                    expected=f"bytes {start}-{min(end,size-1)}/{size}"
                    if r.headers["Content-Range"]!=expected or (self.size is not None and size!=self.size):
                        raise ValueError("Inconsistent FITS HTTP range")
                    data=r.content
                if len(data)!=min(end,size-1)-start+1:
                    raise ValueError("Incomplete FITS range")
                self.size=size
                self.receipts.append([start,len(data),sha256(data).hexdigest()])
                return data
            except (OSError,ValueError,KeyError):
                if attempt:
                    raise
                time.sleep(.2)

    def seekable(self):
        return True

    def readable(self):
        return True

    def seek(self,offset,whence=0):
        self.pos=(0 if whence==0 else self.pos if whence==1 else self.size)+offset
        return self.pos

    def tell(self):
        return self.pos

    def read(self,n=-1):
        stop=self.size if n<0 else min(self.size,self.pos+n)
        if stop<=self.pos:
            return b""
        # Header probes, fibermap and a single row in each wavelength/flux
        # HDU. Overlapping reads reuse exact blocks. No entire HEALPix FITS
        # is retained simply to examine one matched object.
        block=131072
        blocks=[]
        while self.pos<stop:
            start=self.pos//block*block
            if start not in self.parts:
                self.parts[start]=self.fetch(start,min(start+block-1,self.size-1))
            end=min(start+block,stop)
            blocks.append(self.parts[start][self.pos-start:end-start])
            self.pos=end
        return b"".join(blocks)

    def close(self):
        self.session.close()
        super().close()


def fetch_spectrum(m,cache):
    """Cache a single validated spectrum. All public responses get a digest."""
    key=m.get("sparcl_id") or coadd_url(m)+"?TARGETID="+m["specid"]
    file=cache/(sha256(key.encode()).hexdigest()+".npz")
    receipt=file.with_suffix(".json")
    if file.exists() and receipt.exists():
        return file,json.loads(receipt.read_text())
    log=dict(key=key,attempts=[])
    try:
        if m.get("sparcl_id"):
            if not m["data_release"].startswith(("DESI-","BOSS-","SDSS-")):
                return None,dict(log,status="dataset_not_supported")
            params=dict(format="json",units="False",include="sparcl_id,specid,wavelength,flux,ivar,mask",
                        dataset_list=m["data_release"])
            for a in range(2):
                try:
                    raw=requests.post(SPARCL,json=[key],params=params,timeout=(5,26))
                    raw.raise_for_status()
                    data=raw.json()
                    rec=next(x for x in data if isinstance(x,dict) and str(x.get("sparcl_id"))==key)
                    w,f,iv=[np.asarray(rec[k],dtype=float) for k in ("wavelength","flux","ivar")]
                    mask=np.asarray(rec["mask"],dtype=float)
                    log["sha256"]=sha256(raw.content).hexdigest()
                    log["attempts"].append(dict(status="read",url=SPARCL,bytes=len(raw.content)))
                    break
                except (OSError,ValueError,TypeError,KeyError,StopIteration) as exc:
                    log["attempts"].append(dict(status="failed",url=SPARCL,error=str(exc)[:240]))
                    if a==1:
                        raise ValueError("SPARCL spectrum unavailable on both attempts") from exc
                    time.sleep(.5)
        else:
            url=coadd_url(m)
            mirror=url.replace("https://data.desi.lbl.gov/","https://portal.nersc.gov/cfs/desi/")
            for endpoint in (url,mirror):
                try:
                    with RangeFile(endpoint) as rf:
                        with fits.open(rf,memmap=False) as hd:
                            ids=hd["FIBERMAP"].data["TARGETID"]
                            row=np.flatnonzero(ids==int(m["specid"]))
                            if len(row)!=1:
                                raise ValueError("DESI FIBERMAP did not identify exactly one target")
                            idx=int(row[0])
                            parts=[]
                            for band in "BRZ":
                                wave=hd[band+"_WAVELENGTH"].data
                                vals,ivar,mask=[hd[band+"_"+name].section[idx,:] for name in ("FLUX","IVAR","MASK")]
                                parts.append((np.array(wave),np.array(vals),np.array(ivar),np.array(mask)))
                        log["attempts"].append(dict(status="read",url=endpoint,ranges=rf.receipts,targetid=m["specid"]))
                    break
                except (OSError,KeyError,ValueError,IndexError) as exc:
                    log["attempts"].append(dict(url=endpoint,status="failed",error=str(exc)[:200]))
                    if endpoint==mirror:
                        raise
            wc,fc,ic,mc=[np.concatenate([p[k] for p in parts]) for k in range(4)]
            # DESI's camera-arm overlap shares a calibrated 0.8 Å sampling
            # grid. Combine only identical wavelengths with independent IVAR.
            w,inverse=np.unique(wc.round(5),return_inverse=True)
            usable=np.isfinite(fc)&np.isfinite(ic)&(ic>0)&(mc==0)
            iv=np.bincount(inverse,weights=np.where(usable,ic,0.))
            f=np.divide(np.bincount(inverse,weights=np.where(usable,fc*ic,0.)),iv,
                        out=np.zeros_like(w),where=iv>0)
            mask=~(iv>0)
        if not (len(w)==len(f)==len(iv)==len(mask) and len(w)>100 and np.all(np.diff(w)>0)):
            raise ValueError("invalid external spectral wavelength grid")
        e=np.divide(1.,np.sqrt(np.maximum(iv,0)),out=np.full_like(iv,np.nan,dtype=float),
                    where=(iv>0)&(mask==0))
        np.savez_compressed(file,w=w,f=f,e=e)
        log.update(status="measured",array_sha256=sha256(file.read_bytes()).hexdigest())
        receipt.write_text(dump(log))
        return file,log
    except (OSError,ValueError,TypeError,KeyError,IndexError) as exc:
        log["attempts"].append(dict(status="unavailable",error=str(exc)[:300]))
        return None,dict(log,status="spectrum_unavailable",error=str(exc)[:300])


def batched_matches(found):
    positions=sorted({(r["peak"]["ra"],r["peak"]["dec"]) for r in found})
    tasks=[(positions[i:i+BATCH],s) for i in range(0,len(positions),BATCH)
           for s in ("SPARCL","CDS_DESI_DR1")]
    bypoint=defaultdict(dict)
    receipts=[]
    def run(task):
        block,service=task
        return block,service,*query_catalog(block,service)
    with ThreadPoolExecutor(max_workers=4) as pool:
        for block,service,objects,r in pool.map(run,tasks):
            batch=len(receipts)
            receipts.append(r)
            for coord in block:
                outcome=dict(status="unavailable" if objects is None else "queried",batch=batch,matches=[])
                for m in objects or []:
                    sep=float(separation(*coord,number(m["ra"]),number(m["dec"])))
                    if np.isfinite(sep) and sep<=RADIUS:
                        outcome["matches"].append(dict(m,sep_arcsec=sep,catalog=service))
                bypoint[coord][service]=outcome
    return bypoint,receipts


def trace_spectra(found,bypoint,cache):
    # Query sky positions once, retrieve each spectrum at most once. Rerun
    # forced lines at each independent wavelength hypothesis, never copying
    # diagnostics at a displaced center from a neighbor's search.
    cache.mkdir(parents=True,exist_ok=True)
    rows=[m for p in bypoint.values() for s in p.values() for m in s["matches"]]
    unique={}
    for m in rows:
        key=m.get("sparcl_id") or coadd_url(m)+"?TARGETID="+m["specid"]
        m["spectrum_key"]=key
        unique[key]=m
    def obtain(item):
        key,m=item
        return key,fetch_spectrum(m,cache)
    ids=[(k,m) for k,m in unique.items() if m.get("sparcl_id")]
    with ThreadPoolExecutor(max_workers=4) as pool:
        spectra=dict(pool.map(obtain,ids))
    # SPARCL can hold several coadds. For a uniquely represented DESI DR1
    # target reuse the retrieved arrays, retaining both catalogue receipts.
    dr1=defaultdict(list)
    for k,m in ids:
        if m["data_release"].upper()=="DESI-DR1" and spectra[k][0] is not None:
            dr1[m["specid"]].append(k)
    mirror=[(k,m) for k,m in unique.items() if not m.get("sparcl_id")]
    counts=Counter(m["specid"] for _,m in mirror)
    todo=[]
    for k,m in mirror:
        siblings=dr1.get(m["specid"],[])
        if len(siblings)==counts[m["specid"]]==1:
            path,log=spectra[siblings[0]]
            spectra[k]=(path,dict(status="measured",same_unique_target_as=siblings[0],
                                  array_sha256=log["array_sha256"]))
        else:
            todo.append((k,m))
    with ThreadPoolExecutor(max_workers=3) as pool:
        spectra.update(pool.map(obtain,todo))
    fits_cache={}
    for r in found:
        p=r["peak"]
        svc=bypoint[(p["ra"],p["dec"])]
        ext=r["external"]=dict(radius_arcsec=RADIUS,services={s:{k:v for k,v in o.items() if k!="matches"}
                                                             for s,o in svc.items()},matches=[])
        for s in svc.values():
            for m in s["matches"]:
                path,log=spectra[m["spectrum_key"]]
                obj=dict(m,status=log["status"],spectral_receipt_key=m["spectrum_key"])
                if path:
                    key=(str(path),p["wave"],p["channel"])
                    if key not in fits_cache:
                        try:
                            with np.load(path) as a:
                                fits_cache[key]=measure_arrays(a["w"],a["f"],a["e"],p["wave"],p["channel"])
                        except (ValueError,IndexError,TypeError,np.linalg.LinAlgError) as exc:
                            fits_cache[key]=dict(error=str(exc)[:240])
                    diag=fits_cache[key]
                    if "error" in diag:
                        obj.update(status="fit_failed",reason=diag["error"])
                    else:
                        obj["diagnostics"]=diag
                        ls=diag.get("lines",{})
                        # Measured at two independently specified wavelengths
                        # in one external spectrum. These are descriptive fixed
                        # centroids, not maximized over the line-fit grid.
                        obj["same_spectrum_uv_pair"]=(number(ls.get("lya",{}).get("snr"))>=5 and
                              number(ls.get("heii",{}).get("snr"))>=4)
                ext["matches"].append(obj)
        ext["status"]=("matched" if ext["matches"] else "no_counterpart"
                       if all(q["status"]=="queried" for q in svc.values()) else "partial_or_unavailable")
        classify_matches(r)
    return {k:receipt for k,(_,receipt) in spectra.items()}


def classify_matches(r):
    p,ext=r["peak"],r["external"]
    z=row_z(p)
    anchored=[m for m in ext["matches"] if number(m["redshift_warning"])==0 and
              number(m["redshift"])>=0 and abs(number(m["redshift"])-z)<.006]
    r["independent_redshift_agrees"]=bool(anchored)
    r["external_fixed_uv_pair"]=any(m.get("same_spectrum_uv_pair") for m in ext["matches"])
    r["first_VAC_agrees"]="independent DESI visual redshift" in r["first"].get("passed",[])
    r["measured_spectra"]=sum("diagnostics" in m for m in ext["matches"])
    r["selection_reason"]=("external_redshift" if anchored else "both_uv_lines_in_external_spectrum" if r["external_fixed_uv_pair"]
                           else "external_spectrum" if r["measured_spectra"]
                           else "independent_two_line_repeat" if r["same_visit_two_line_repeats"]
                           else "earlier_DESI_visual_redshift" if r["first_VAC_agrees"]
                           else "catalog_match_only" if ext["matches"] else ext["status"])


def focus(indir,out,cache):
    indir,out,cache=map(Path,(indir,out,cache))
    out.mkdir(parents=True,exist_ok=True)
    found,total=select(indir)
    bypoint,tap_receipts=batched_matches(found)
    spec_receipts=trace_spectra(found,bypoint,cache/"focus-spectra")
    return write_catalog(indir,out,found,total,tap_receipts,spec_receipts)


def write_catalog(indir,out,found,total,tap_receipts,spec_receipts):
    fields=("row_index","channel","stream","wave","detectid","source_id","shotid","ra","dec")
    header=fields+("first_pass_line","wave_hex","trial_redshift","group_id","first_pass_external","confirmed_repeat_shots",
                   "SPARCL_query","CDS_query","matches","measured_spectral_entries","external_redshift_agrees",
                   "prior_VAC_agrees","external_fixed_uv_pair","followup_status","spectroscopic_ids","eligibility","static_DESI","ASU_DESI")
    with gzip.open(out/"catalog.csv.gz","wt",newline="") as f, gzip.open(out/"followup.jsonl.gz","wt") as g:
        writer=csv.DictWriter(f,header)
        writer.writeheader()
        for r in found:
            p=r["peak"]
            ext=r["external"]
            status=r["selection_reason"]
            total["eligible"]+=1
            total["discovery_"+status]+=1
            m=ext["matches"]
            writer.writerow(dict({k:p.get(k) for k in fields},first_pass_line=r["first_pass_line"],wave_hex=r["key"][-1],trial_redshift=row_z(p),
                            group_id=r["group"],first_pass_external=r["first_pass_external"],
                            confirmed_repeat_shots=dump(r["same_visit_two_line_repeats"]),
                            SPARCL_query=ext["services"]["SPARCL"]["status"],
                            CDS_query=ext["services"]["CDS_DESI_DR1"]["status"],matches=len(m),
                            measured_spectral_entries=r["measured_spectra"],external_redshift_agrees=r["independent_redshift_agrees"],
                            prior_VAC_agrees=r["first_VAC_agrees"],external_fixed_uv_pair=r["external_fixed_uv_pair"],followup_status=status,
                            spectroscopic_ids=dump([x["spectrum_key"] for x in m]),eligibility=dump(r["criteria"]),
                           static_DESI=ext["services"].get("DESI_DR1_static",{}).get("status"),
                           ASU_DESI=ext["services"].get("CDS_ASU_DESI_DR1",{}).get("status")))
            if m or r["same_visit_two_line_repeats"] or r["first_VAC_agrees"]:
                g.write(dump(r)+"\n")
                total["followup_records"]+=1
    total["physical_groups"]=len({r["group"] for r in found})
    total["distinct_query_positions"]=len({(r["peak"]["ra"],r["peak"]["dec"]) for r in found})
    summary=dict(counts=dict(total),input_sha256=sha256((indir/"vetted.jsonl.gz").read_bytes()).hexdigest(),
        scope="All unresolved HeII lines with Lyα ≥5σ, observed HeII/Lyα ≥0.1, FWHM <900 km/s, "
        "no reported ≥3σ metal component; unknown coverage retained. All unresolved spectral "
        "or imaging EW0 lower bounds >240 Å. Genuine source streams only. 1.25 arcsec cones.",
        rule="Discovery table status is a crossmatch result. Only complete adjudications in vetted.jsonl.gz "
        "can revise the first-pass outcome. Stable per-release key: row_index, channel, stream, wave.hex().",
        catalogs=["NOIRLab SPARCL (DESI, BOSS, SDSS)","CDS V/161/zcatdr1 (DESI DR1 Redrock)"])
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    with gzip.open(out/"receipts.json.gz","wt") as f:
        f.write(dump(dict(catalog=tap_receipts,spectra=spec_receipts)))
    return summary


def static_retry(indir,out,release,cache,vet):
    """Retry failed positions with independent public DESI files; leave successful requests intact."""
    from desi_static import search as static_search
    from evidence_summary import measurement_id
    from acquisition import Spectra,prepare
    keep=("catalog.csv.gz","followup.jsonl.gz","summary.json","receipts.json.gz","vetted.jsonl.gz",
          "vet-summary.json","vet-inputs.json","inputs-"+release+".json")
    archive=out/"first-attempt"
    if not archive.exists():
        with tempfile.TemporaryDirectory(dir=out) as temp:
            snapshot=Path(temp)/"first-attempt"
            snapshot.mkdir()
            for name in keep:
                if (out/name).is_file():
                    shutil.copy2(out/name,snapshot/name)
            for name in keep[:6]:
                if not (snapshot/name).is_file():
                    raise ValueError("Missing first-attempt evidence: "+name)
            snapshot.rename(archive)
    original=json.loads((archive/"summary.json").read_text())
    digest=sha256((indir/"vetted.jsonl.gz").read_bytes()).hexdigest()
    if digest!=original["input_sha256"]:
        raise ValueError("Saved focus products belong to a different HETDEX ledger")
    done=out/"static-summary.json"
    if done.exists():
        meta=json.loads(done.read_text())
        if meta["input_sha256"]!=digest:
            raise ValueError("Previous static retry belongs to different inputs")
        if not meta["position_status"].get("partial_or_unavailable") and (not vet or meta["requested_vet"]):
            print("Previously completed static retry",meta["position_status"],flush=True)
            return
    base,total=select(indir)
    with gzip.open(archive/"catalog.csv.gz","rt") as stream:
        prior_csv=list(csv.DictReader(stream))
    if len(base)!=len(prior_csv):
        raise ValueError("Previously eligible catalog differs from frozen selection")
    old=list(lines(archive/"followup.jsonl.gz"))
    restore_legacy_line_numbers(old,indir/"vetted.jsonl.gz")
    carried={r["first_pass_line"]:r for r in old}
    old_vet=list(lines(archive/"vetted.jsonl.gz"))
    # Historical full vetting used this stable order. Validate the unchanged
    # measured primary plus partner identity, group, and exact original decision.
    sortkey=lambda r:(r["first_pass_external"],not r["independent_redshift_agrees"],
                      not bool(r["same_visit_two_line_repeats"]),r["group"])
    if len(old)!=len(old_vet):
        raise ValueError("First-attempt follow-up and adjudication counts disagree")
    old_order=iter(sorted(old,key=sortkey))
    decisions={}
    for item in old_vet:
        parent=next(old_order)
        if (measurement_id(item["peak"])!=measurement_id(parent["peak"]) or
            item["group"]!=parent["group"] or
            item.get("first_pass_line",parent["first_pass_line"])!=parent["first_pass_line"]):
            raise ValueError("Historical focus decision does not map uniquely to the original source row")
        item["first_pass_line"]=parent["first_pass_line"]
        decisions[item["first_pass_line"]]=item
    challenged=[]
    for r,csvrow in zip(base,prior_csv,strict=True):
        p=r["peak"]
        if (row_key(p)!=(int(csvrow["row_index"]),csvrow["channel"],csvrow["stream"],csvrow["wave_hex"]) or
            (p["ra"],p["dec"])!=(float(csvrow["ra"]),float(csvrow["dec"]))):
            raise ValueError("First-attempt catalog order disagrees with frozen source measurements")
        if r["first_pass_line"] in carried:
            prior=carried[r["first_pass_line"]]
            r["external"]=prior["external"]
        else:
            r["external"]=dict(radius_arcsec=RADIUS,matches=[],status=csvrow["followup_status"],
                       services={s:dict(status=csvrow[c]) for s,c in (("SPARCL","SPARCL_query"),("CDS_DESI_DR1","CDS_query"))})
        if (csvrow["SPARCL_query"]!="queried" or csvrow["CDS_query"]!="queried" or
            int(csvrow["matches"]) and not int(csvrow["measured_spectral_entries"])):
            challenged.append(r)
        classify_matches(r)
    coords=[(r["peak"]["ra"],r["peak"]["dec"]) for r in challenged]
    cone,log=static_search(coords)
    bypoint={p:{"static_retry":dict(status=c["status"],matches=c["matches"])} for p,c in cone.items()}
    fresh=[dict(r) for r in challenged]
    new_spectra=trace_spectra(fresh,bypoint,cache/"focus-spectra")
    revet=[]
    for r,new in zip(challenged,fresh,strict=True):
        p=r["peak"]
        c=cone[(p["ra"],p["dec"])]
        ext=r["external"]
        ext["services"].update(DESI_DR1_static=dict(status=c["static_status"],files=c["files"],
                    healpixels=c["expected_HEALPix"]),CDS_ASU_DESI_DR1=dict(status=c["ASU_status"]))
        prior={m["spectrum_key"]:m for m in ext["matches"]}
        changed=False
        for m in new["external"]["matches"]:
            key=m["spectrum_key"]
            if key not in prior or (not prior[key].get("diagnostics") and m.get("diagnostics")):
                prior[key]=m
                changed=True
        ext["matches"]=list(prior.values())
        ext["status"]="matched" if prior else c["status"]
        if not prior and c["status"]=="partial_or_unavailable" and ext["services"]["CDS_DESI_DR1"]["status"]=="queried":
            ext["status"]="no_desi_dr1_counterpart"
        classify_matches(r)
        if changed:
            revet.append(r)
    status=Counter(c["status"] for c in cone.values())
    with tempfile.TemporaryDirectory(dir=out) as temp:
        stage=Path(temp)
        with gzip.open(archive/"receipts.json.gz","rt") as f:
            original_receipts=json.load(f)
        summary=write_catalog(indir,stage,base,total,original_receipts["catalog"],
                              {**original_receipts["spectra"],**new_spectra})
        summary["catalogs"].extend(["DESI DR1 native iron Redrock files","CDS single-cone V/161/zcatdr1 ASU"])
        summary["first_attempt"]="first-attempt/"
        summary["static_retry"]="static-summary.json"
        (stage/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
        with gzip.open(stage/"static-receipts.json.gz","wt") as f:
            f.write(dump(log))
        if vet and revet:
            sub=stage/"new"
            sub.mkdir()
            with gzip.open(sub/"followup.jsonl.gz","wt") as f:
                for r in revet:
                    f.write(dump(r)+"\n")
            path=prepare(release,cache,sub)
            with Spectra(path,release) as cat:
                adjudicate(cat,sub,cache,indir)
            for rr in lines(sub/"vetted.jsonl.gz"):
                decisions[rr["first_pass_line"]]=rr
            shutil.copy2(sub/"vet-inputs.json",stage/"static-vet-inputs.json")
            shutil.copy2(sub/("inputs-"+release+".json"),stage/("static-inputs-"+release+".json"))
        with gzip.open(stage/"vetted.jsonl.gz","wt") as f:
            for n in sorted(decisions):
                f.write(dump(decisions[n])+"\n")
        counts=Counter(r["vet"]["disposition"] for r in decisions.values())
        (stage/"vet-summary.json").write_text(json.dumps(dict(rows=len(decisions),counts=dict(counts),
                new_counterpart_adjudications=len(revet) if vet else 0),indent=2)+"\n")
        stat=dict(input_sha256=digest,requested_vet=vet,eligible_rows=len(base),retried_rows=len(challenged),
                  searched_positions=len(cone),position_status=dict(status),new_matches_for_vetting=len(revet),
                  new_adjudications=len(revet) if vet else 0,coadd_receipts=len(new_spectra),
                  archive_file_status=dict(Counter(r["status"] for r in log["redrock"])),
                  static_position_status=dict(Counter(x["static_status"] for x in cone.values())),
                  single_cone_status=dict(Counter(x["ASU_status"] for x in cone.values())),
                  scope="Static empty cones apply only to DESI DR1 released products. Unavailable SDSS/BOSS "
                  "remains untested. No published HEALPix products, successful target-cone searches and "
                  "request failures are distinct in catalog.csv.gz.")
        for f in stage.iterdir():
            if f.is_file():
                f.replace(out/f.name)
        done.write_text(json.dumps(stat,indent=2)+"\n")
    print("DESI static fallback",json.dumps(stat,sort_keys=True),flush=True)


def restore_legacy_line_numbers(records,original):
    """Match complete frozen peak and original decision; duplicate identical records keep ordinal order."""
    pending=defaultdict(deque)
    def token(p,v):
        raw=json.dumps([p,v],sort_keys=True,separators=(",",":"),allow_nan=False)
        return sha256(raw.encode()).digest()
    for r in records:
        if not r.get("first_pass_line"):
            pending[token(r["peak"],r["first"])].append(r)
    if not pending:
        return
    for n,base in enumerate(lines(original),1):
        key=token(base["peak"],base["vet"])
        if key in pending:
            pending[key].popleft()["first_pass_line"]=n
            if not pending[key]:
                del pending[key]
            if not pending:
                return
    raise ValueError("Saved follow-up rows did not match the original ledger exactly")


def adjudicate(cat,out,cache_root,indir=None):
    """Run the full original veto chain on every matched or repeat-supported row, with no rank cap."""
    from vetting import Ancillary, Cache, enrich, external_check, tap_batch, validate_local
    cache=Cache(cache_root)
    aux=Ancillary(cache)
    counts=Counter()
    def prepare(r):
        p=dict(r["peak"])
        c=cat.metadata(int(p["row_index"]))
        p={**c,**p,"catalog":c}
        p.update(enrich(p,cat))
        pair=p.get("pair_lya") or {}
        if number(pair.get("snr"))>number((p.get("lya_velocity_best") or {}).get("snr")):
            p["lya_velocity_best"]=dict(pair,air_wave=pair["wave"],offset_km_s=p.get("pair_velocity_kms"))
        v=validate_local(p,cat,aux)
        if aux.errors.get("LAN"):
            v["LAN_catalog_unavailable"]=True
        return p,v
    records=list(lines(out/"followup.jsonl.gz"))
    restore_legacy_line_numbers(records,(indir or out.parent)/"vetted.jsonl.gz")
    # The frozen first pass checked separation and uniqueness of positional
    # repeat visits as well as shared source IDs. Restore those same row pairs.
    cat.confirmed_repeat_rows={}
    for item in records:
        idx=item["peak"]["row_index"]
        cat.confirmed_repeat_rows[idx]=tuple(m["row_index"] for m in item["first"].get("repeats",[])
                                             if m.get("row_index")!=idx and m.get("shotid")!=item["peak"].get("shotid"))
    # More independent support first, but the loop has no numerical shortlist limit.
    records.sort(key=lambda r:(r["first_pass_external"],not r["independent_redshift_agrees"],
                               not bool(r["same_visit_two_line_repeats"]),r["group"]))
    try:
        with gzip.open(out/"vetted.jsonl.gz","wt") as f:
            for i in range(0,len(records),12):
                block=records[i:i+12]
                ready=[]
                for r in block:
                    p,v=prepare(r)
                    ready.append((p,v,r))
                selected=[(p,v,r) for p,v,r in ready if v["disposition"]=="unresolved"]
                if selected:
                    points=[row[0] for row in selected]
                    with ThreadPoolExecutor(max_workers=3) as pool:
                        batches=list(pool.map(lambda s:tap_batch(s,points),("NED","SIMBAD","Legacy")))
                    for j,(p,v,r) in enumerate(selected):
                        for s,answer in zip(("NED","SIMBAD","Legacy"),batches):
                            v[s]=answer[j]
                        supplies=dict(p.get("rescue_audit") or {})
                        supplies["SPARCL"]=r["external"]
                        for k in ("HSC","Legacy_forced","cube"):
                            prior=r["first"].get(k) or {}
                            if prior.get("status")=="measured" and k not in supplies:
                                supplies[k]=prior
                        p["rescue_audit"]=supplies
                    def final_check(t):
                        p,v,r=t
                        return external_check(p,v,cache)
                    with ThreadPoolExecutor(max_workers=3) as pool:
                        for (_,v,_),answer in zip(selected,pool.map(final_check,selected)):
                            v.update(answer)
                for p,v,r in ready:
                    counts[v["disposition"]]+=1
                    p.pop("catalog",None)
                    p.pop("rescue_audit",None)
                    f.write(dump(dict(peak=p,vet=v,key=r["key"],first_pass_line=r["first_pass_line"],group=r["group"],
                            first_pass_disposition=r["first"]["disposition"],
                            first_pass_external=r["first_pass_external"]))+"\n")
                    f.flush()
    finally:
        aux.close()
    (out/"vet-summary.json").write_text(json.dumps(dict(rows=len(records),counts=dict(counts),
           note="Every selected ledger hypothesis restored from the original release and rechecked."),indent=2)+"\n")
    (out/"vet-inputs.json").write_text(dump(cache.inputs)+"\n")


ID_REASONS={
    "heii":{"independent line identification unavailable"},
    "ew":{"external redshift unavailable for low-z and AGN discrimination",
          "independent Ly-alpha identification required for imaging-led EW"},
    "blob":{"independent redshift unavailable"}}


def replay_pair_decision(p,v,restore_visits):
    """Revisit only identification-dependent gates in already fully vetted real rows."""
    from vetting import AGN_GAP, AGN_PASSED, PAIR_EXTERNAL, PAIR_REPEAT, fixed_external_pairs
    from vetting import ordinary_agn_disfavored, stack_pair_shots
    ch=p["channel"]
    if (ch not in ID_REASONS or v.get("pair_evidence_version")==2 or "SPARCL" not in v or
        v["disposition"] not in ("candidate","unresolved")):
        return None
    old=set(v.get("reasons",[]))
    passed=list(v.get("passed",[]))
    pairs=fixed_external_pairs(p,v)
    stack=p.get("stream")=="repeat_stack"
    primary=(PAIR_EXTERNAL in passed or "independent instrument confirms line" in passed or
             "independent HETDEX shot" in passed)
    anchor=bool({"independent DESI visual redshift","independent survey redshift agrees"} & set(passed))
    relevant=ID_REASONS[ch] | {"foreground interpretation unresolved"}
    reviewed=(ch=="ew" and p.get("stream")=="imaging_rescue")
    lowwave=ch in ("ew","blob") and p["wave"]<3714 and "HETDEX Ly-alpha redshift agrees" in passed
    benefited_from_stack=stack and not anchor and not pairs and (
        PAIR_REPEAT in passed or (lowwave and not (old & relevant) and primary))
    if not (pairs and old & relevant or benefited_from_stack):
        return None
    v=dict(v)
    reasons=set(old)
    if pairs:
        passed=list(dict.fromkeys(passed+[PAIR_EXTERNAL]))
        v["external_fixed_uv_pairs"]=pairs
        reasons.difference_update(relevant)
    if benefited_from_stack:
        shots=stack_pair_shots(p,v)
        if len(shots)<2:
            recovered=restore_visits(p,v)
            if recovered is not None:
                v["secondary_visit_measurements"]=recovered
                shots=stack_pair_shots(p,v)
        v["stack_pair_shots"]=shots
        if len(shots)<2:
            passed=[text for text in passed if text!=PAIR_REPEAT]
            if ch=="heii":
                reasons.add("independent line identification unavailable")
            elif ch=="ew":
                reasons.add("external redshift unavailable for low-z and AGN discrimination")
                if reviewed:
                    reasons.add("independent Ly-alpha identification required for imaging-led EW")
            else:
                reasons.add("independent redshift unavailable")
            if v.get("catalog_conflict") or v.get("low_z_families"):
                reasons.add("foreground interpretation unresolved")
    # In the original pipeline the AGN gate ran only after all other gates
    # cleared. A newly identified source must still pass it. Preserve all
    # existing independent imaging, morphology, foreground and novelty vetoes.
    if not reasons:
        if ordinary_agn_disfavored(p,v,row_z(p)):
            passed=list(dict.fromkeys(passed+[AGN_PASSED]))
        else:
            reasons.add(AGN_GAP)
    outcome="unresolved" if reasons else "candidate"
    if reasons==old and outcome==v["disposition"]:
        return None
    return dict(v,disposition=outcome,reasons=sorted(reasons),passed=passed,pair_evidence_version=2)


def cached_secondary_visits(indir,cache,stack,receipts):
    """Optional missing reference-visit flux; never fetch absent files."""
    from acquisition import Spectra,checksums
    from vetting import repeat_line
    release=indir.name
    opened=[]
    def restore(p,v):
        if not opened:
            cat=None
            manifest=indir/("inputs-"+release+".json")
            if cache is not None and manifest.is_file():
                known=json.loads(manifest.read_text())
                files=[(cache/r["file"],r["sha256"]) for r in known]
                if len(files)==2 and all(q.is_file() for q,_ in files):
                    for path,expected in files:
                        if checksums(path)["sha256"]!=expected:
                            raise ValueError("Cached HETDEX input disagrees with original release SHA-256")
                    spec=next(q for q,_ in files if "_spec_" in q.name)
                    cat=stack.enter_context(Spectra(spec,release))
                    receipts["verified_cached_files"]={q.name:h for q,h in files}
            opened.append(cat)
        if opened[0] is None:
            receipts["source_visits_unavailable"]+=1
            return None
        if p["channel"]=="heii":
            other=p.get("lya_velocity_best") or {}
        else:
            fam=p.get("lya_family") or {}
            other=max((x for k,x in fam.items() if k in ("heii","oiii_uv1","oiii_uv2","ciii1","ciii2")),
                      key=lambda x:number(x.get("snr")),default={})
        rows=p.get("repeat_rows") or []
        if not rows or number(other.get("snr"))<5 or number(other.get("air_wave"))<=0:
            receipts["precise_secondary_or_visits_unavailable"]+=1
            return None
        receipts["remeasured_stack_secondary"]+=1
        return repeat_line(opened[0],int(p["row_index"]),{**p,"wave":other["air_wave"]},rows=rows)
    return restore


def replay_saved(indir,out,cache=None):
    """Generate an immutable, original-line-addressed overlay; performs no network requests."""
    from evidence_summary import adjudications
    out.mkdir(parents=True,exist_ok=True)
    paths=["vetted.jsonl.gz","focus/vetted.jsonl.gz","repeat-blobs/vetted.jsonl.gz"]
    def digest(path):
        with path.open("rb") as f:
            from hashlib import file_digest
            return file_digest(f,"sha256").hexdigest()
    hashes={name:digest(indir/name) for name in paths if (indir/name).is_file()}
    count=Counter()
    receipts=Counter()
    target=out/"saved-pair-decisions.jsonl.gz"
    tmp=target.with_suffix(".tmp")
    with ExitStack() as stack:
        restore=cached_secondary_visits(indir,cache,stack,receipts)
        with gzip.open(tmp,"wt") as f:
            for ordinal,(line,record,source,_) in enumerate(adjudications(indir,include_pair_replay=False),1):
                p,v=record["peak"],record["vet"]
                count["effective_rows_examined"]+=1
                current=replay_pair_decision(p,v,restore)
                if current is None:
                    continue
                count["revised_rows"]+=1
                count[v["disposition"]+"_to_"+current["disposition"]]+=1
                f.write(dump(dict(peak=p,vet=current,key=row_key(p),first_pass_line=ordinal,
                             previous_ledger=source,previous_ledger_line=line,previous_disposition=v["disposition"]))+"\n")
    tmp.replace(target)
    summary=dict(input_sha256=hashes["vetted.jsonl.gz"],source_ledgers=hashes,
                 counts=dict(count),cached_visit_receipts=dict(receipts),
                 method="Fixed Ly-alpha ≥5σ and HeII ≥4σ in one separate external spectrum, and two distinct "
                 "individual primary/secondary visits for stack-only pair confirmation. All unrelated vetoes "
                 "are preserved. Missing member spectra leave the repeat criterion unresolved.")
    (out/"saved-pair-summary.json").write_text(json.dumps(summary,indent=2)+"\n")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("--cache",type=Path,default=Path(os.environ.get("RUNNER_TEMP","/tmp"))/"hetdex")
    p.add_argument("--release",choices=("sc1","sc2"))
    p.add_argument("--vet",action="store_true")
    p.add_argument("--reuse-matches",action="store_true",help="Use checksummed focus outputs already present.")
    p.add_argument("--static-retry",action="store_true",help="Retry unavailable first-pass spectroscopic cones in the public DESI DR1 HEALPix products.")
    p.add_argument("--replay-saved",action="store_true",help="Recheck earlier adjudications from saved evidence without network access.")
    args=p.parse_args()
    if args.static_retry:
        static_retry(args.input,args.output,args.release or args.input.name,args.cache,args.vet)
        return
    if args.replay_saved:
        replay_saved(args.input,args.output,args.cache)
        return
    if not args.reuse_matches:
        current=focus(args.input,args.output,args.cache)
    else:
        current=json.loads((args.output/"summary.json").read_text())
        digest=sha256((args.input/"vetted.jsonl.gz").read_bytes()).hexdigest()
        if current["input_sha256"]!=digest:
            p.error("Previous focus inputs do not match the provided vetted ledger")
    if args.vet:
        release=args.release or args.input.name
        if release not in ("sc1","sc2"):
            p.error("A source release is needed to restore the HETDEX spectra")
        if current["counts"].get("followup_records",0)==0:
            with gzip.open(args.output/"vetted.jsonl.gz","wt"):
                pass
            (args.output/"vet-summary.json").write_text('{"rows":0,"counts":{}}\n')
            return
        from acquisition import Spectra,prepare
        path=prepare(release,args.cache,args.output)
        with Spectra(path,release) as cat:
            adjudicate(cat,args.output,args.cache,args.input)


if __name__=="__main__":
    main()
