#!/usr/bin/env python3
"""Independent official multiwavelength audit of the FROZEN eRASS1 main screen.

No optical non-detection is inferred from the absence of an association row.
Uses new June-2026 Salvato et al. 2025 DR1 LS10/Gaia DR3/CatWISE2020 products.
Network and all computation are performed on the single GitHub Actions job.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import sys
import tarfile
import tempfile
import traceback
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits

DIR = Path(__file__).resolve().parent
OUT = DIR / "output"
ROOT = "https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/SalvatoM_DR1/"
RELEASES = {
  "LS10": "Salvato_etal2025_DR1_LS10.colfits.tgz",
  "GDR3": "Salvato_etal2025_DR1_GDR3.colfits.tgz",
  "CW2020": "Salvato_etal2025_DR1_CW2020.colfits.tgz",
}
# The thresholds are fixed before this script's first execution, not optimized to the observed associations.
PROTOCOL = {
  "release": "Salvato et al. 2025 (A&A 704 A344) official eRASS1 optical/IR counterparts, distributed 2026-06-09",
  "surveys": list(RELEASES), "secure_NWAY_p_any":0.8, "tentative_NWAY_p_any":0.1,
  "reference_main":"eRASS1 Main.v1.2 full-screen outputs from run.py",
  "join":"exact X-ray UID (not 8-arcsec sky proximity); choose max NWAY p_any × p_i per catalogue",
  "limit_definition":"5-sigma LS10 g-band point-source depth from PSFDEPTH_G only where positive and in footprint",
  "fx_fopt_proxy":"Maccacaro log10(f_X/f_V)=log10(f_X)+0.4*V+5.37, substituted g for V for context ONLY (uncalibrated)",
  "do_not_conclude":"no association != no object; LS10 depths not automatically available at absent-counterpart positions",
}
RECEIPT = {"began_utc":datetime.now(timezone.utc).isoformat(),
           "protocol":PROTOCOL,"downloads":[],"status":"started"}
AUDIT = {"sources":0, "top60":0,"tables":{}, "scope":"frozen main-screen sources, direct official counterpart joins"}
def sha(path):
    dig=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b""): dig.update(chunk)
    return dig.hexdigest()

def save(name,data):
    (OUT/name).write_text(json.dumps(data,sort_keys=True,indent=2,allow_nan=False)+"\n")

def download(session,name,url,dest):
    rec={"name":name,"url":url,"time_utc":datetime.now(timezone.utc).isoformat()}
    for retry in (1,2):
        try:
            with session.get(url,stream=True,timeout=(30,180)) as r:
                rec.update({"http":r.status_code,"resolved_url":r.url,
                            "last_modified":r.headers.get("Last-Modified"),
                            "etag":r.headers.get("ETag"),"attempt":retry,
                            "content_type":r.headers.get("Content-Type")})
                r.raise_for_status()
                with open(dest,"wb") as f:
                    for chunk in r.iter_content(8*1024*1024):
                        if chunk: f.write(chunk)
            rec.update({"archive_bytes":dest.stat().st_size,"archive_sha256":sha(dest)})
            RECEIPT["downloads"].append(rec)
            return rec
        except Exception as e:
            rec["error"]=repr(e)
    RECEIPT["downloads"].append(rec)
    raise RuntimeError(name+" unavailable: "+str(rec))

def extract(tgz,tmp,name,rec):
    with tarfile.open(tgz,"r:*") as z:
        fs=[m for m in z.getmembers() if m.isfile() and m.name.lower().endswith((".colfits",".fits"))]
        if len(fs)!=1: raise RuntimeError(name+" expected one COLFITS member, got "+repr([(m.name,m.size) for m in fs]))
        m=fs[0]
        path=tmp/(name+".colfits")
        with z.extractfile(m) as src,open(path,"wb") as out:
            for chunk in iter(lambda:src.read(8*1024*1024),b""):out.write(chunk)
    rec.update({"fits_member":m.name,"fits_bytes":path.stat().st_size,"fits_sha256":sha(path)})
    return path

def column(tbl,name):
    if name not in tbl.names:return None
    ar=np.asarray(tbl[name])
    if ar.ndim>1 and ar.shape[0]==1: ar=ar[0]
    return ar.reshape(-1)
def at(cols,k,i):
    a=cols.get(k)
    if a is None: return None
    v=a[i]
    if isinstance(v,np.generic):v=v.item()
    if isinstance(v,bytes):return v.decode("utf-8","replace").strip()
    if isinstance(v,str):return v.strip()
    if isinstance(v,(float,np.floating)):
        return float(v) if math.isfinite(float(v)) and abs(float(v))<1e18 else None
    if isinstance(v,(bool,np.bool_)):return bool(v)
    if isinstance(v,(int,np.integer)):return int(v)
    return None

def read_csv(name):
    with (OUT/name).open(newline="") as f:return list(csv.DictReader(f))
def write_csv(name, rows, fields):
    with (OUT/name).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore")
        w.writeheader()
        for row in rows:w.writerow({k:row.get(k,"") for k in fields})

EXTRAS={
 "LS10":["NWAY_ncat","NWAY_p_any","NWAY_p_i","NWAY_match_flag","NWAY_Separation_LS10_ERO",
         "LS10_OBJID","LS10_FLUX_G","LS10_FLUX_R","LS10_FLUX_W1","LS10_FLUX_W2",
         "PSFDEPTH_G","InAnyLS10","InAllLS10","has_Salvato18_AGN_colors",
         "has_Salvato22_AGN_colors","class_gal_exgal","coronal","p_stellar",
         "Gaia_moving_5sigma","zphot","zspec_compilation","exgal","bestCTP",
         "with_Quaia","with_GaiaQSO","with_R90"],
 "GDR3":["NWAY_ncat","NWAY_p_any","NWAY_p_i","NWAY_match_flag",
         "NWAY_Separation_GDR3_ERO","GDR3_source_id_1","parallax","pmra","pmdec",
         "phot_g_mean_mag","parallax_over_error","GDR3_source_id_1a"],
 "CW2020":["NWAY_ncat","NWAY_p_any","NWAY_p_i","NWAY_match_flag",
           "NWAY_Separation_CW2020_ERO","CW2020_source_id",
           "CW2020_w1mpro","CW2020_w2mpro"],
}
def run():
    rows=read_csv("screened_main.csv")
    top=read_csv("candidates.csv")
    uidset={int(r["uid"]) for r in rows}
    assert len(uidset)==len(rows),"duplicate source UID in screened main"
    AUDIT["sources"]=len(rows);AUDIT["top60"]=len(top)
    by_uid={int(r["uid"]):{"uid":r["uid"],"iauname":r["iauname"],
                             "ra":r["ra"],"dec":r["dec"],
                             "fX_0p2_2p3":r["flux_0p2_2p3"],"counts":r["ml_counts"],
                             "g_equals_V_required_for_1e4":r["illustrative_g_equals_V_needed_for_ratio_1e4"]}
            for r in rows}
    success=0
    with requests.Session() as session, tempfile.TemporaryDirectory(prefix="erass1-counterparts-") as d:
        session.headers.update({"User-Agent":"Open-Space-XDINS-official-counterpart-audit/1.0"})
        tmp=Path(d)
        for name,filename in RELEASES.items():
            url=ROOT+filename
            try:
                tgz=tmp/(name+".tgz")
                rec=download(session,name,url,tgz)
                colfits=extract(tgz,tmp,name,rec)
                with fits.open(colfits,memmap=True,lazy_load_hdus=True) as f:
                    tbl=f[1].data
                    keys=["UID",*EXTRAS[name]]
                    cols={key:column(tbl,key) for key in keys}
                    if cols["UID"] is None:raise RuntimeError("UID missing in "+name)
                    n=len(cols["UID"])
                    if n < 900000:raise RuntimeError(name+" has implausible row count "+str(n))
                    rec.update({"rows_or_vector_length":n,"hdu_rows":len(tbl),
                                "hdu_name":f[1].name,"missing_requested_fields":[k for k,v in cols.items() if v is None]})
                    selected=np.flatnonzero(np.isin(cols["UID"], list(uidset)))
                    per_source={}
                    for i in selected:
                        u=int(cols["UID"][i])
                        candidate={k:at(cols,k,i) for k in EXTRAS[name]}
                        is_real=candidate.get("NWAY_ncat")==2
                        p=candidate.get("NWAY_p_any")
                        q=candidate.get("NWAY_p_i")
                        candidate["has_catalogue_counterpart"]=bool(is_real)
                        candidate["secure_counterpart"]=bool(is_real and p is not None and p>=PROTOCOL["secure_NWAY_p_any"])
                        score=(float(p) if isinstance(p,float) else 0)*(
                            float(q) if isinstance(q,float) else 0)
                        # Every UID has at most one selected "best" match in this reduced table;
                        # retain count of all alternatives and choose best probabilities.
                        entry=per_source.setdefault(u,{"n_joined_rows":0,"best_score":-1,"best":None})
                        entry["n_joined_rows"]+=1
                        if (score>entry["best_score"] or entry["best"] is None):
                            entry.update({"best_score":score,"best":candidate})
                    for u,entry in per_source.items():
                        by_uid[u][name]=entry["best"]
                        by_uid[u][name+"_rows"]=entry["n_joined_rows"]
                rec["status"]="ok"
                success+=1
                AUDIT["tables"][name]={"status":"ok","rows_in_release":rec["rows_or_vector_length"],
                     "matched_xray_UIDs":len(per_source),
                     "sources_with_association":sum(x["best"]["has_catalogue_counterpart"] for x in per_source.values()),
                     "sources_with_secure_association":sum(x["best"]["secure_counterpart"] for x in per_source.values()),
                     "missing_fields":rec["missing_requested_fields"]}
                colfits.unlink(missing_ok=True);tgz.unlink(missing_ok=True)
            except Exception as e:
                for r in RECEIPT["downloads"]:
                    if r["name"]==name:r["status"]="failed_parse_or_download";r["error"]=repr(e)
                AUDIT["tables"][name]={"status":"failed","reason":repr(e)}
    follow=[]
    for u,d in by_uid.items():
        l=d.get("LS10",{});g=d.get("GDR3",{});c=d.get("CW2020",{})
        d["ls10_associated"]=l.get("has_catalogue_counterpart",None)
        d["ls10_secure"]=l.get("secure_counterpart",None)
        d["gaia_associated"]=g.get("has_catalogue_counterpart",None)
        d["gaia_secure"]=g.get("secure_counterpart",None)
        d["catwise_associated"]=c.get("has_catalogue_counterpart",None)
        d["catwise_secure"]=c.get("secure_counterpart",None)
        d["LS10_p_any"]=l.get("NWAY_p_any")
        d["GDR3_p_any"]=g.get("NWAY_p_any")
        d["CW2020_p_any"]=c.get("NWAY_p_any")
        d["ls10_agn_colors"]=bool(l.get("has_Salvato18_AGN_colors") or l.get("has_Salvato22_AGN_colors")) if l else None
        d["ls10_coronal"]=l.get("coronal")
        d["ls10_zphot"]=l.get("zphot")
        d["ls10_exgal"]=l.get("exgal")
        d["ls10_stellar_prob"]=l.get("p_stellar")
        d["ls10_g_mag"]=None
        d["ls10_g_5sigma_limit"]=None
        d["g_proxy_fx_fopt_lower_limit"]=None
        if l:
            fg=l.get("LS10_FLUX_G")
            if fg is not None and fg>0:d["ls10_g_mag"]=round(22.5-2.5*math.log10(fg),3)
            depth=l.get("PSFDEPTH_G")
            if depth is not None and depth>0:
                glimit=22.5-2.5*math.log10(5/math.sqrt(depth))
                d["ls10_g_5sigma_limit"]=round(glimit,3)
                if not l.get("has_catalogue_counterpart",False):
                    # g vs V conversion is not certified: only a descriptive, unvalidated proxy.
                    value=float(d["fX_0p2_2p3"])
                    if value>0:d["g_proxy_fx_fopt_lower_limit"]=round(10**(math.log10(value)+.4*glimit+5.37),1)
        strong=[s for s in RELEASES if d.get({"LS10":"ls10_secure","GDR3":"gaia_secure","CW2020":"catwise_secure"}[s])]
        d["secure_official_identifications"]=";".join(strong)
        d["status"]="OFFICIAL_SECURE_COUNTERPART" if strong else (
            "NO_SECURE_ASSOCIATION_LIMITS_UNVERIFIED" if success==len(RELEASES) else "INCOMPLETE_OFFICIAL_ASSOCIATIONS")
        follow.append(d)
    follow.sort(key=lambda r:r["uid"])
    cols=[
       "uid","iauname","ra","dec","fX_0p2_2p3","counts",
       "ls10_associated","ls10_secure","LS10_p_any","ls10_g_mag","ls10_g_5sigma_limit",
       "ls10_agn_colors","ls10_coronal","ls10_zphot","ls10_exgal","ls10_stellar_prob",
       "gaia_associated","gaia_secure","GDR3_p_any","catwise_associated","catwise_secure","CW2020_p_any",
       "secure_official_identifications","g_equals_V_required_for_1e4",
       "g_proxy_fx_fopt_lower_limit","status"]
    write_csv("official_counterparts.csv",follow,cols)
    topuid={x["uid"] for x in top}
    write_csv("official_top60.csv",[x for x in follow if x["uid"] in topuid],cols)
    counts={"source_rows":len(follow),
            "top60_count":len(topuid),
            "any_secure_counterpart_in_full_screen":sum(bool(x["secure_official_identifications"]) for x in follow),
            "no_secure_counterpart_in_full_screen":sum(not bool(x["secure_official_identifications"]) for x in follow),
            "top60_secure_counterpart":sum(bool(x["secure_official_identifications"]) for x in follow if x["uid"] in topuid),
            "top60_no_secure_counterpart":sum(not bool(x["secure_official_identifications"]) for x in follow if x["uid"] in topuid),
            "top60_ls10_5sigma_g_depth_available":sum(x["ls10_g_5sigma_limit"] is not None for x in follow if x["uid"] in topuid),
            "all_sources_validated_xdins":0}
    AUDIT["counts"]=counts
    AUDIT["assertions"]=[
        "Association tables are source-matching probabilities, NOT 5-sigma optical non-detections.",
        "A missing LS10 join row or absent counterpart cannot establish that Legacy Surveys covered the exact position.",
        "Even with local g depth, Maccacaro V-based flux ratio with g substituted is not physically calibrated.",
        "The screened 1,330 are not a complete unbiased soft-X-ray XDINS sample; native spectra/timing and pixel checks remain absent.",
    ]
    AUDIT["pending"]=[
        "Pan-STARRS1, DES, unWISE and 2MASS independent detection+limiting-depth tests at 1,330 positions",
        "Image-level LS DR10 g/r/z non-detection upper limits and coverage for no-counterpart sky positions",
        "Native XMM, Chandra, Swift spectra, light curves and pulsation searches",
        "Verified Kurpas+ other 2024-26 prior-object matching and latest ATNF radio pulsar catalogue",
        "Independent completeness injection and variability rejection",
    ]
    AUDIT["release_tables_completed"]=success
    RECEIPT["status"]="all_official_tables_checked" if success==len(RELEASES) else ("partial_official_tables" if success else "failed_all_official_tables")
    save("official_followup.json",AUDIT)
    check=json.loads((OUT/"checks.json").read_text())
    check["counts"].update({"official_"+k:v for k,v in counts.items()})
    check["official_followup"]={"status":RECEIPT["status"],"tables":AUDIT["tables"]}
    for name in RELEASES:
        if AUDIT["tables"].get(name,{}).get("status")=="ok":
            check["completed"].append("official eRASS1 "+name+" release by X-ray UID (not standalone per-image depth)")
        else:check["pending"].append("official eRASS1 "+name+" counterpart download/parse")
    save("checks.json",check)
    if not success:raise RuntimeError("all official tables failed; not a null result")

if __name__=="__main__":
    try:run()
    except Exception as e:
        RECEIPT["error"]=repr(e)
        RECEIPT["traceback"]=traceback.format_exc()[-3800:]
        if RECEIPT["status"]=="started":RECEIPT["status"]="failed"
    finally:
        RECEIPT["finished_utc"]=datetime.now(timezone.utc).isoformat()
        save("official_receipts.json",RECEIPT)
        print(json.dumps({"status":RECEIPT["status"],"tables":AUDIT.get("tables"),"counts":AUDIT.get("counts"),"error":RECEIPT.get("error")}))
    if RECEIPT["status"] in ("started","failed","failed_all_official_tables"):sys.exit(2)
