"""Uncapped second-pass spectroscopy of eligible, previously unresolved HETDEX hypotheses."""
import argparse
from collections import Counter, defaultdict
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
    for row in lines(indir/"vetted.jsonl.gz"):
        p,v=row["peak"],row["vet"]
        counters["ledger_rows"]+=1
        criteria=eligible(p,v)
        if criteria is None:
            continue
        counters[p["channel"]+"_eligible"]+=1
        found.append(dict(peak=p,first=v,criteria=criteria,
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
            with RangeFile(url) as rf:
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
                log["attempts"].append(dict(status="read",url=url,ranges=rf.receipts,targetid=m["specid"]))
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
    return {k:receipt for k,(_,receipt) in spectra.items()}


def focus(indir,out,cache):
    indir,out,cache=map(Path,(indir,out,cache))
    out.mkdir(parents=True,exist_ok=True)
    found,total=select(indir)
    bypoint,tap_receipts=batched_matches(found)
    spec_receipts=trace_spectra(found,bypoint,cache/"focus-spectra")
    fields=("row_index","channel","stream","wave","detectid","source_id","shotid","ra","dec")
    header=fields+("wave_hex","trial_redshift","group_id","first_pass_external","confirmed_repeat_shots",
                   "SPARCL_query","CDS_query","matches","measured_spectral_entries","external_redshift_agrees",
                   "prior_VAC_agrees","external_fixed_uv_pair","followup_status","spectroscopic_ids","eligibility")
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
            writer.writerow(dict({k:p.get(k) for k in fields},wave_hex=r["key"][-1],trial_redshift=row_z(p),
                            group_id=r["group"],first_pass_external=r["first_pass_external"],
                            confirmed_repeat_shots=dump(r["same_visit_two_line_repeats"]),
                            SPARCL_query=ext["services"]["SPARCL"]["status"],
                            CDS_query=ext["services"]["CDS_DESI_DR1"]["status"],matches=len(m),
                            measured_spectral_entries=r["measured_spectra"],external_redshift_agrees=r["independent_redshift_agrees"],
                            prior_VAC_agrees=r["first_VAC_agrees"],external_fixed_uv_pair=r["external_fixed_uv_pair"],followup_status=status,
                            spectroscopic_ids=dump([x["spectrum_key"] for x in m]),eligibility=dump(r["criteria"])))
            if m or r["same_visit_two_line_repeats"] or r["first_VAC_agrees"]:
                g.write(dump(r)+"\n")
                total["followup_records"]+=1
    total["physical_groups"]=len({r["group"] for r in found})
    total["distinct_query_positions"]=len(bypoint)
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


def adjudicate(cat,out,cache_root):
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
                    f.write(dump(dict(peak=p,vet=v,key=r["key"],group=r["group"],
                            first_pass_disposition=r["first"]["disposition"],
                            first_pass_external=r["first_pass_external"]))+"\n")
                    f.flush()
    finally:
        aux.close()
    (out/"vet-summary.json").write_text(json.dumps(dict(rows=len(records),counts=dict(counts),
           note="Every selected ledger hypothesis restored from the original release and rechecked."),indent=2)+"\n")
    (out/"vet-inputs.json").write_text(dump(cache.inputs)+"\n")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("--cache",type=Path,default=Path(os.environ.get("RUNNER_TEMP","/tmp"))/"hetdex")
    p.add_argument("--release",choices=("sc1","sc2"))
    p.add_argument("--vet",action="store_true")
    p.add_argument("--reuse-matches",action="store_true",help="Use checksummed focus outputs already present.")
    args=p.parse_args()
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
            adjudicate(cat,args.output,args.cache)


if __name__=="__main__":
    main()
