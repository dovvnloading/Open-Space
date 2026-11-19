#!/usr/bin/env python3
"""Bounded recovery of unscreened frozen DR14 pn sources from HEASARC PPS.

First retry improved coordinate-verified PNX000 merged pn products. If unavailable,
search PN per-exposure SRCTSR8 products, but mark them PARTIAL and do not
promote them to complete observation screens. Every FITS match is astrometric.
No absent PPS product is a QPE non-detection.
"""
import datetime as dt
import hashlib
import html
import json
import re
import time
from pathlib import Path
from urllib.parse import urljoin, urlparse, unquote

import numpy as np
import requests
import control_validation as cv
from heasarc_batch import (OUT, BASE, HOST, CUTOFF, SESSION, manifest, save,
  get_file, net_merged, flare_background, background_peak_comparison,
  background_only_finder, poisson_injection, scan_observation, utc)

def pn_listing(obsid, rec):
    url=HOST+obsid+"/PPS/"
    for attempt in range(1,4):
        try:
            with SESSION.get(url,timeout=(12,30)) as resp:
                resp.raise_for_status()
                text=resp.text
            names=[]
            for href in re.findall(r'href\s*=\s*["\x27]([^"\x27]+)',text,flags=re.I):
                parsed=urlparse(urljoin(url,html.unescape(href)))
                if parsed.netloc.lower()!="heasarc.gsfc.nasa.gov":
                    continue
                if not parsed.path.startswith("/FTP/xmm/data/rev0/"+obsid+"/PPS/"):
                    continue
                name=unquote(parsed.path.rsplit("/",1)[-1])
                if re.fullmatch(r"P"+obsid+r"PN[A-Z0-9]{4}SRCTSR8[0-9A-F]{3}\.FTZ",name,re.I):
                    names.append(name)
            names=sorted(set(names))
            rec.setdefault("directory_retries",[]).append({
                "attempt":attempt,"http_status":200,"n_pn_src_tsr":len(names)})
            return names
        except Exception as exc:
            rec.setdefault("directory_retries",[]).append({
                "attempt":attempt,"error":repr(exc)[:300]})
            if attempt<3:time.sleep(attempt)
    return []

def per_exposure_fallback(obsid, needed, record, limit=60):
    rows={r["detid"]:r for r in needed}
    files=[n for n in pn_listing(obsid,record) if "PNX000" not in n.upper()]
    if not files:
        record["per_exposure_status"]="no_pn_per_exposure_SRCTSR8"
        return []
    expected={f"{int(r['pps_srcnum']):03X}" for r in needed if r.get("pps_srcnum","").isdigit()}
    preferred=[f for f in files if f[-7:-4].upper() in expected]
    nonpreferred=[f for f in files if f not in preferred]
    picked=(preferred+nonpreferred)[:limit]
    record["per_exposure_available"]=len(files)
    record["per_exposure_products_selected"]=len(picked)
    record["per_exposure_uninspected"]=len(files)-len(picked)
    # Source-masked field particle backgrounds are separate from source
    # extraction BACK1-3 but not equivalent to soft-band event photometry.
    bgdata=None
    bgfiles=record.get("directory_receipts",[])
    if bgfiles:
        # Use one source-masked field background in exact same observation.
        # scan_observation() already wrote a SHA for the earlier download.
        fb=record.get("field_particle_background",{})
        if fb.get("filename"):
            raw=get_file(obsid,fb["filename"],record["file_receipts"],attempts=2)
            if raw is not None:
                val=flare_background(raw,fb["filename"])
                if isinstance(val,tuple): _,bgdata=val
    rescued={}
    for fn in picked:
        if len(rescued)==len(needed):
            break
        raw=get_file(obsid,fn,record["file_receipts"],attempts=2)
        if raw is None: continue
        try:
            pos,curve,ref,hdr=net_merged(raw)
        except Exception as exc:
            record.setdefault("per_exposure_parse_errors",[]).append({
                "file":fn,"error":repr(exc)[:240]})
            continue
        distances=sorted((cv.sep_arcsec(pos,(float(r["ra"]),float(r["dec"]))),r["detid"])
                         for r in needed)
        nearest,detid=distances[0]
        if nearest>15 or (len(distances)>1 and distances[1][0]-nearest<2):
            continue
        if detid in rescued:continue
        finder=cv.find_eruptions(curve,ref,(.5,24.))
        result={
            "obs_id":obsid,"detid":detid,"srcid":rows[detid]["srcid"],
            "status":"position_verified_pn_single_exposure_net_soft_screened_PARTIAL",
            "coverage":"one pn exposure product, not all observation GTIs",
            "source_sep_arcsec":round(nearest,3),
            "pn_product":fn,"pn_fits_sha256":hashlib.sha256(raw).hexdigest(),
            "source_fits_header":hdr,
            "good_300s_bins":int(np.count_nonzero(curve["good"])),
            "finder":finder,
            "unvetted_repeated_peak_trigger":bool(finder.get("range_match") and finder.get("n_peaks",0)>=2),
            "independent_field_particle_background":background_peak_comparison(finder,bgdata),
            "offsource_scaled_background_control":background_only_finder(curve),
            "poisson_count_domain_injection":poisson_injection(curve,detid),
            "pileup_SAS_epatplot":"NOT PERFORMED",
            "soft_offsource_event_region":"NOT PERFORMED",
            "validated_new_QPE":False,
        }
        rescued[detid]=result
    record["per_exposure_status"]="PARTIAL_PRODUCT_SCREEN" if rescued else "no_coordinate_verified_pn_exposure"
    return list(rescued.values())

def known_screened():
    manifest_rows,freeze=manifest()
    confirmed=set()
    for path in sorted(OUT.glob("heasarc_batch_cohort_*.json")):
        j=json.loads(path.read_text())
        if j.get("cohort_manifest_sha256")!=freeze["final_csv_sha256"]:
            raise ValueError("Cohort receipt hash changed: "+path.name)
        for o in j.get("observations",[]):
            for s in o.get("source_results",[]):
                if s.get("status")=="position_verified_net_pn_soft_screened":
                    confirmed.add((o["obs_id"],s["detid"]))
    for path in sorted(OUT.glob("heasarc_rescue_*.json")):
        j=json.loads(path.read_text())
        if j.get("manifest_sha256")!=freeze["final_csv_sha256"]:
            raise ValueError("Rescue receipt cohort hash changed: "+path.name)
        for o in j.get("observations",[]):
            for s in o.get("source_results",[])+o.get("partial_source_results",[]):
                if s.get("status")=="position_verified_net_pn_soft_screened":
                    confirmed.add((o["obs_id"],s["detid"]))
    return manifest_rows,freeze,confirmed

def main():
    cfg=json.loads((BASE/"control_batch.json").read_text())
    assert cfg["stage"]=="heasarc_rescue"
    batch=str(cfg["batch_id"])
    assert re.fullmatch(r"[a-zA-Z0-9_]+",batch)
    requested=cfg["target_obsids"]
    assert isinstance(requested,list) and 1<=len(requested)<=180
    assert len(requested)==len(set(requested))
    cohort,fr,old=known_screened()
    assert all(oid in cohort for oid in requested)
    rec={"batch_id":batch,"started_utc":utc(),"stage":"frozen_missing_pn_source_retries",
         "source":"NASA HEASARC mirror PPS, exact pn product receipts",
         "manifest_sha256":fr["final_csv_sha256"],
         "finder_threshold_sha256":cv.checksum(BASE/"control_thresholds.json"),
         "target_obsids":requested,"observations":[],"status":"RUNNING",
         "partial_exposures_not_counted_as_complete":True}
    name=f"heasarc_rescue_{batch}.json";t0=time.monotonic()
    for oid in requested:
        if time.monotonic()-t0>CUTOFF:
            rec["status"]="TIME_BUDGET_PAUSE_NOT_COMPLETE";break
        needed=[r for r in cohort[oid] if (oid,r["detid"]) not in old]
        if not needed:
            row={"obs_id":oid,"status":"already_screened","source_results":[],"partial_source_results":[]}
        else:
            row=scan_observation(oid,needed)
            completed={s["detid"] for s in row["source_results"]
                       if s.get("status")=="position_verified_net_pn_soft_screened"}
            uncompleted=[r for r in needed if r["detid"] not in completed]
            row["partial_source_results"]=per_exposure_fallback(oid,uncompleted,row) if uncompleted else []
            for s in row["source_results"]:
                if s.get("status")=="position_verified_net_pn_soft_screened":
                    old.add((oid,s["detid"]))
        rec["observations"].append(row)
        save(name,rec)
        print("HEASARC_RESCUE",batch,oid,row["status"],
              "merged",sum(x.get("status")=="position_verified_net_pn_soft_screened" for x in row["source_results"]),
              "partial",len(row["partial_source_results"]),flush=True)
    if rec["status"]=="RUNNING":rec["status"]="COMPLETE_WITH_REVIEW_LIMITS"
    rec["elapsed_s"]=round(time.monotonic()-t0,2)
    rec["finished_utc"]=utc()
    rec["limitations"]=[
      "Missing sources are re-tried with coordinate verified merged pn products (full screen) first",
      "If no merged product exists, a single per-exposure pn time series is PARTIAL, not a negative or complete observation",
      "No calibrated event-list source/exposure union, SAS pile-up, independent soft event-region extraction or blind archive FPR",
      "Source-specific BACK1-3 and source-masked field particle light curves cannot alone confirm a QPE",
      "Crossmatched galaxy host sample is SDSS/Gaia provisional, not a complete stellar/4LAC/Legacy-survey host census",
    ]
    save(name,rec)
    from heasarc_postrescue import main as referee_main
    referee_main()
    # Trigger follow-ups are optional and independently receipt-logged.
    # Retrieval failures must never erase the completed pn source batch.
    from heasarc_trigger_vet import main as vet_main
    from heasarc_trigger_curves import main as curve_main
    for stage, function in (("catalogs", vet_main), ("curves", curve_main)):
        try:
            function()
        except Exception as exc:
            save(f"heasarc_rescue_auxiliary_failure_{stage}.json",
                 {"status": "INCOMPLETE", "stage": stage,
                  "error": repr(exc)[:600],
                  "frozen_QPE_validation": False})
    print("HEASARC_RESCUE_SUMMARY",batch,rec["status"],len(rec["observations"]),flush=True)

if __name__=="__main__":
    main()
