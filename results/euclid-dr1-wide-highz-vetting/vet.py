#!/usr/bin/env python3
"""Conservative empirical colour-model and bounded external vetting of full-scan survivors.

No external service failure is a nondetection. A Q-versus-dwarf model is always
conditional on its training classes and explicitly declared prior; no galaxies,
artifacts, or footprint-incompleteness are included.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import time
from xml.etree import ElementTree as ET

import numpy as np
import requests

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
MODEL_PRIOR = 0.001
TOP_EXTERNAL = 24
ARCSEC = 1.0/3600.0

def utc():
    return datetime.now(timezone.utc).isoformat()

def save(path, x):
    path.write_text(json.dumps(x, indent=2, allow_nan=False, sort_keys=True)+"\n")

def good(v):
    try:
        f = float(v)
        return f if math.isfinite(f) else None
    except (ValueError, TypeError):
        return None

def vector_from_control(hit):
    if not hit:
        return None
    y = good(hit.get("Y_minus_J_AB"))
    h = good(hit.get("J_minus_H_AB"))
    if y is None or h is None or not (-5 <= y <= 8 and -5 <= h <= 5):
        return None
    return (y, h)

def model_training(controls):
    q, d, identifiers = [], [], {"quasar": [], "dwarf": []}
    for x in controls:
        v = vector_from_control(x.get("nearest"))
        if v is not None:
            group = "quasar" if x["class"] == "quasar" else "dwarf"
            (q if group == "quasar" else d).append(v)
            identifiers[group].append(x["name"])
    return np.array(q, dtype=float).reshape(-1,2), np.array(d, dtype=float).reshape(-1,2), identifiers

def mean_likelihood(x, points, bandwidth=.32):
    # Explicit fixed-bandwidth bivariate empirical kernel density;
    # much larger than most catalogue high-SNR random colour errors.
    if not len(points):
        return 0.0
    squared = np.sum(((points-np.asarray(x))/bandwidth)**2,axis=1)
    return float(np.mean(np.exp(-.5*squared))) / (2*math.pi*bandwidth**2)

def p_q_given_q_or_d(x, q, d, prior=MODEL_PRIOR):
    lq,ld = mean_likelihood(x,q), mean_likelihood(x,d)
    den = prior*lq + (1-prior)*ld
    return float(prior*lq/den) if den>0 else None

def diagnostics(q,d):
    result = {"matched_clean_colours_quasar":len(q),
              "matched_clean_colours_dwarf":len(d),
              "frozen_colour_bandwidth_mag":.32,
              "conditional_priors": [.001,.01],
              "model_validated_for_population": False,
              "training_set_selection":"Yang et al. 2026 known quasars and spectroscopic dwarf contaminants matched within 3 arcsec",
              "omitted_classes":["red galaxies","artifacts","blends","other stellar classes"],
              "method":"empirical 2-D (Y-J,J-H) fixed-bandwidth Gaussian KDE; conditional Bayesian reweighting only"}
    if len(q)<5 or len(d)<5:
        result["model_status"]="INSUFFICIENT_MATCHED_LABELED_CONTROLS"
        return result
    tests=[]
    for target,same,other,is_q in ((p,q,d,True) for p in []):
        pass
    for group, x, other, truth in (("quasar",q,d,True),("dwarf",d,q,False)):
        for i, point in enumerate(x):
            leave=np.delete(x,i,axis=0)
            pred = p_q_given_q_or_d(point, leave, other, prior=.001) if truth else p_q_given_q_or_d(point,other,leave,prior=.001)
            tests.append({"group":group,"p":pred,"truth":int(truth)})
    result["crossvalidated_n"]=len(tests)
    result["crossvalidated_mean_p_for_quasars"]=round(sum(v["p"] for v in tests if v["truth"])/len(q),5)
    result["crossvalidated_mean_p_for_dwarfs"]=round(sum(v["p"] for v in tests if not v["truth"])/len(d),5)
    result["crossvalidated_brier_score"]=round(sum((v["p"]-v["truth"])**2 for v in tests)/len(tests),6)
    result["model_status"]="UNVALIDATED_CONDITIONAL_TWO_CLASS_MODEL"
    result["model_caution"]="Published controls are discovery-selected rather than a population-representative held-out sample; numerical posteriors are not empirical purity or discovery significance."
    return result

def parse_csv(response):
    if response.status_code != 200:
        raise RuntimeError("HTTP_"+str(response.status_code)+" "+response.text[:150])
    txt=response.text.lstrip("\ufeff")
    if txt.lstrip().startswith("<"):
        try:
            root=ET.fromstring(response.content)
            details=[{"name":x.get("name"),"value":x.get("value"),"message":(x.text or "")[:200]}
                     for x in root.findall(".//{*}INFO")]
            raise RuntimeError("TAP_VOTABLE_RESPONSE_"+str(details)[:320])
        except ET.ParseError:
            raise RuntimeError("SERVICE_XML_UNPARSEABLE_"+txt[:160].replace("\n"," "))
    lines=[line for line in txt.splitlines()
           if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        return []
    rows=list(csv.DictReader(lines))
    if rows and "ra" not in rows[0] and "s_ra" not in rows[0] and "main_id" not in rows[0]:
        raise RuntimeError("UNEXPECTED_TABLE_SCHEMA:"+str(list(rows[0])[:8]))
    return rows

def sphere_distance(ra1,dec1,ra2,dec2):
    a,b=math.radians(dec1), math.radians(dec2)
    da=math.radians(ra1-ra2)
    val=math.sin((a-b)/2)**2+math.cos(a)*math.cos(b)*math.sin(da/2)**2
    return math.degrees(2*math.asin(math.sqrt(max(0,min(1,val)))))*3600

def tap_request(endpoint, sql, label, fields, r, timeout=16):
    record={"service":label,"endpoint":endpoint,"adql":sql,"utc":utc(),"status":"QUERY_ERROR",
            "catalog_match":None,"nearest_arcsec":None}
    try:
        response=requests.get(endpoint,params={"REQUEST":"doQuery","LANG":"ADQL",
                              "FORMAT":"csv","QUERY":sql},timeout=(8,timeout),
                              headers={"User-Agent":"Open-Space-Euclid-DR1F-vetting/1.0"})
        record["http_status"]=response.status_code
        record["response_sha256"]=hashlib.sha256(response.content).hexdigest()
        record["response_bytes"]=len(response.content)
        matches=parse_csv(response)
        record["returned_rows"]=len(matches)
        record["status"]="CATALOG_ROWS" if matches else "NO_CATALOG_ROWS_NOT_A_FLUX_LIMIT"
        chosen=None
        for m in matches:
            ra=good(m.get(fields[0]));dec=good(m.get(fields[1]))
            if ra is None or dec is None:
                continue
            sep=sphere_distance(float(r["ra"]),float(r["dec"]),ra,dec)
            if sep<=3 and (chosen is None or sep<chosen[0]):
                chosen=(sep,m)
        if chosen:
            record["catalog_match"]=chosen[1]
            record["nearest_arcsec"]=round(chosen[0],4)
    except Exception as ex:
        record["error"]=type(ex).__name__+": "+str(ex)[:220]
    return record

def query_irsa(r, table, fields):
    ra,dec=float(r["ra"]),float(r["dec"])
    sql=("SELECT TOP 10 "+",".join(fields)+" FROM "+table+
         " WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',"
         f"{ra:.8f},{dec:.8f},0.001))")
    return tap_request("https://irsa.ipac.caltech.edu/TAP/sync",sql,table,("ra","dec"),r)

def query_gaia(r):
    ra,dec=float(r["ra"]),float(r["dec"])
    sql=("SELECT TOP 10 source_id,ra,dec,pmra,pmra_error,pmdec,pmdec_error,"
         "parallax,parallax_error FROM gaiadr3.gaia_source "
         "WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
         f"CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.00075))")
    return tap_request("https://gea.esac.esa.int/tap-server/tap/sync",sql,
                       "Gaia_DR3",("ra","dec"),r)

def query_simbad(r):
    ra,dec=float(r["ra"]),float(r["dec"])
    sql=("SELECT TOP 10 b.main_id,b.ra,b.dec,b.otype FROM basic AS b "
         "WHERE 1=CONTAINS(POINT('ICRS',b.ra,b.dec),"
         f"CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.001))")
    # CDS announced a Nov 2026 Strasbourg outage. Query an explicit
    # independently hosted mirror, preserving its ingestion lag, then fail over.
    attempts=[]
    for endpoint,asof in (
        ("https://simbad.cfa.harvard.edu/simbad/sim-tap/sync","2026-10-14"),
        ("https://simbad.idia.ac.za/simbad/sim-tap/sync","UNVERIFIED")):
        item=tap_request(endpoint,sql,"SIMBAD_MIRROR_literature_not_exhaustive",
                         ("ra","dec"),r,timeout=18)
        item["mirror_ingest_date"]=asof
        if item["status"]!="QUERY_ERROR":
            item["prior_failed_mirrors"]=attempts
            return item
        attempts.append({"endpoint":endpoint,"error":item.get("error"),
                         "http_status":item.get("http_status")})
    return {"service":"SIMBAD_MIRRORS","status":"QUERY_ERROR",
            "mirror_attempts":attempts,"catalog_match":None,
            "limitations":"CDS main site outage and mirror ingestion lag; unknown literature is not a nondetection"}

def query_spectra(r,source):
    ra,dec=float(r["ra"]),float(r["dec"])
    if source=="DESI_DR1":
        # NOIRLab TAP's ADQL layer does not support CONTAINS/CIRCLE in this
        # database. Use an indexable bounding box and postfilter to 3 arcsec.
        half_dec=.00084
        half_ra=half_dec/max(.1,math.cos(math.radians(dec)))
        lo,hi=ra-half_ra,ra+half_ra
        if lo<0:
            ra_condition=f"(ra>={lo+360:.8f} OR ra<={hi:.8f})"
        elif hi>=360:
            ra_condition=f"(ra>={lo:.8f} OR ra<={hi-360:.8f})"
        else:
            ra_condition=f"ra BETWEEN {lo:.8f} AND {hi:.8f}"
        sql=("SELECT TOP 10 id,ra,dec,z FROM desi_dr1.zpix "
             f"WHERE {ra_condition} AND dec BETWEEN {dec-half_dec:.8f} "
             f"AND {dec+half_dec:.8f}")
        return tap_request("https://datalab.noirlab.edu/tap/sync",sql,
                           "DESI_DR1_catalog_redshift_not_spectrum",("ra","dec"),r,timeout=12)
    if source=="ESO":
        sql=("SELECT TOP 10 obs_id,s_ra,s_dec,dataproduct_type FROM ivoa.obscore "
             "WHERE dataproduct_type='spectrum' AND "
             "1=CONTAINS(POINT('ICRS',s_ra,s_dec),"
             f"CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.001))")
        return tap_request("https://archive.eso.org/tap_obs/sync",sql,
                           "ESO_archive_spectrum_metadata",("s_ra","s_dec"),r,timeout=12)
    if source=="SDSS":
        sql=("SELECT TOP 5 s.specObjID,s.z,p.ra,p.dec FROM SpecObj AS s "
             "JOIN PhotoObj AS p ON s.bestObjID=p.objID WHERE "
             f"p.ra BETWEEN {ra-.001:.8f} AND {ra+.001:.8f} "
             f"AND p.dec BETWEEN {dec-.001:.8f} AND {dec+.001:.8f}")
        url="https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
        record={"service":"SDSS_DR18_spectroscopy_metadata","endpoint":url,
                "query":sql,"status":"QUERY_ERROR","utc":utc()}
        try:
            response=requests.get(url,params={"cmd":sql,"format":"csv"},
                                  timeout=(8,12))
            record["http_status"]=response.status_code
            record["response_sha256"]=hashlib.sha256(response.content).hexdigest()
            record["response_bytes"]=len(response.content)
            rows=parse_csv(response)
            record["returned_rows"]=len(rows)
            matches=[]
            for obj in rows:
                qra,qdec=good(obj.get("ra")),good(obj.get("dec"))
                if qra is not None and qdec is not None:
                    sep=sphere_distance(ra,dec,qra,qdec)
                    if sep<=3:
                        matches.append((sep,obj))
            if matches:
                sep,best=min(matches,key=lambda x:x[0])
                record["nearest_arcsec"]=round(sep,4)
                record["catalog_match"]=best
            record["status"]=("MATCHED_SPECTRUM_METADATA_NOT_VALIDATED" if matches
                else "NO_MATCHED_SPECTRUM_METADATA_NOT_A_COMPLETENESS_LIMIT")
        except Exception as ex:
            record["error"]=type(ex).__name__+": "+str(ex)[:220]
        return record
    if source=="MAST":
        url="https://mast.stsci.edu/portal/Mashup/Mashup.asmx/api/v0/invoke"
        payload={"service":"Mast.Caom.Cone","params":{"ra":ra,"dec":dec,"radius":3/3600},
                 "format":"json","pagesize":50}
        record={"service":"MAST_CAOM_all_products","endpoint":url,
                "params":payload,"status":"QUERY_ERROR","utc":utc()}
        try:
            response=requests.post(url,data={"request":json.dumps(payload)},timeout=(8,12))
            record["http_status"]=response.status_code
            record["response_sha256"]=hashlib.sha256(response.content).hexdigest()
            x=response.json()
            if response.status_code!=200 or x.get("status")=="ERROR":
                raise RuntimeError("MAST_"+str(x.get("msg",response.status_code))[:120])
            records=x.get("data") or []
            record["returned_rows"]=len(records)
            record["spectroscopic_products"]=sum(str(t.get("dataproduct_type","")).lower()=="spectrum"
                                                 for t in records)
            record["status"]="METADATA_RETURNED_NEEDS_OBSERVATION_VALIDATION"
        except Exception as ex:
            record["error"]=type(ex).__name__+": "+str(ex)[:220]
        return record
    raise ValueError(source)

def image_metadata(r):
    # Image products only; metadata cannot validate persistence, VIS ghosts or
    # detector edges, which need calibrated per-exposure science cutouts.
    ra,dec=float(r["ra"]),float(r["dec"])
    url="https://irsa.ipac.caltech.edu/SIA"
    record={"service":"Euclid_IRSA_SIA_image_metadata","endpoint":url,
            "status":"QUERY_ERROR","image_level_artifact_status":
            "NOT_EVALUATED_REQUIRES_PIXEL_LEVEL_VIS_AND_NISP_EXPOSURE_REVIEW",
            "utc":utc()}
    try:
        response=requests.get(url,params={
            "COLLECTION":"euclid_dr1_wide_DpdMerBksMosaic",
            "POS":f"CIRCLE {ra:.7f} {dec:.7f} 0.002",
            "FORMAT":"image/fits"},timeout=(8,12))
        record["http_status"]=response.status_code
        record["response_sha256"]=hashlib.sha256(response.content).hexdigest()
        if response.status_code!=200:
            raise RuntimeError("HTTP_"+str(response.status_code))
        record["status"]="IMAGE_METADATA_RESPONSE_NOT_PIXEL_CHECKED"
        record["response_bytes"]=len(response.content)
        record["url_count"]=response.text.count("access_url")
    except Exception as ex:
        record["error"]=type(ex).__name__+": "+str(ex)[:200]
    return record

def motion_and_wise(queries):
    result={"foreground_motion_diagnostic":"NOT_ESTABLISHED",
            "W1_minus_W2_vega_mag":None,
            "wise_color_status":"UNAVAILABLE"}
    wise=queries.get("catwise_2020",{}).get("catalog_match") or {}
    a,b=good(wise.get("w1mpro")),good(wise.get("w2mpro"))
    if a is not None and b is not None:
        result["W1_minus_W2_vega_mag"]=round(a-b,3)
        result["wise_color_status"]="CATWISE_CATALOG_PSF_COLOR_BLENDED_SOURCE_POSSIBLE"
    unwise=queries.get("unwise_2019",{}).get("catalog_match") or {}
    f1,f2=good(unwise.get("flux_1")),good(unwise.get("flux_2"))
    d1,d2=good(unwise.get("dflux_1")),good(unwise.get("dflux_2"))
    result["unwise_W1_over_W2_catalog_flux_ratio"]=(
        round(f1/f2,5) if f1 is not None and f2 is not None and f1>0 and f2>0 else None)
    result["unwise_flux_significances"]= [
        round(a/b,2) if a is not None and b is not None and b>0 else None
        for a,b in ((f1,d1),(f2,d2))]
    result["unwise_color_status"]="FLUX_RATIO_ONLY_NOT_CALIBRATED_VEGA_COLOUR"
    gaia=queries.get("Gaia_DR3",{}).get("catalog_match") or {}
    gaia_sig=[]
    for val,err in (("pmra","pmra_error"),("pmdec","pmdec_error"),
                    ("parallax","parallax_error")):
        v,e=good(gaia.get(val)),good(gaia.get(err))
        if v is not None and e and e>0:
            gaia_sig.append(abs(v/e))
    cw=queries.get("catwise_2020",{}).get("catalog_match") or {}
    cw_sig=[]
    for val,err in (("pmra","sigpmra"),("pmdec","sigpmdec")):
        v,e=good(cw.get(val)),good(cw.get(err))
        if v is not None and e and e>0:
            cw_sig.append(abs(v/e))
    result["gaia_astrometric_component_significances"]= [round(x,2) for x in gaia_sig]
    result["catwise_motion_component_significances"]= [round(x,2) for x in cw_sig]
    if gaia_sig and max(gaia_sig)>5:
        result["foreground_motion_diagnostic"]="REJECT_POSSIBLE_FOREGROUND_GAIA"
    elif cw_sig and max(cw_sig)>5:
        result["foreground_motion_diagnostic"]="REJECT_POSSIBLE_FOREGROUND_CATWISE"
    elif gaia_sig or cw_sig:
        result["foreground_motion_diagnostic"]="NO_SIGNIFICANT_CATALOG_MOTION_NOT_PROOF_EXTRAGALACTIC"
    return result

def vet_one(row, slot):
    queries={}
    for table,fields in (
        ("catwise_2020",["ra","dec","w1mpro","w2mpro","pmra","pmdec","sigpmra","sigpmdec"]),
        ("unwise_2019",["ra","dec","flux_1","dflux_1","flux_2","dflux_2",
                         "fracflux_1","fracflux_2"])):
        queries[table]=query_irsa(row,table,fields)
    queries["Gaia_DR3"]=query_gaia(row)
    queries["SIMBAD"]=query_simbad(row)
    for name in ("DESI_DR1","SDSS","ESO","MAST"):
        queries[name]=query_spectra(row,name)
    if slot < 8:
        queries["Euclid_image_metadata"]=image_metadata(row)
    else:
        queries["Euclid_image_metadata"]={
            "status":"NOT_QUERIED_OUTSIDE_TOP_8",
            "image_level_artifact_status":"NOT_EVALUATED"}
    return {"id":row["id"],"ra":row["ra"],"dec":row["dec"],
            "checked_utc":utc(),"diagnostics":motion_and_wise(queries),
            "queries":queries,
            "image_level_artifact_checked":False,
            "spectroscopic_redshift_confirmed_by_this_study":False,
            "known_quasar_literature_review_complete":False}

def min_training_distance(x, q, d, bandwidth=.32):
    """Distance to nearest labelled control, in fixed KDE bandwidth units."""
    points=np.concatenate((q,d),axis=0)
    return float(np.min(np.linalg.norm((points-np.asarray(x))/bandwidth,axis=1)))


def main():
    summary_file=OUT/"summary.json"
    if not summary_file.exists():
        raise RuntimeError("Aggregate first")
    summary=json.loads(summary_file.read_text())
    if summary["status"]!="SCAN_COMPLETE":
        save(OUT/"vetting_status.json",{"status":"PENDING_FULL_HATS_SCAN",
             "scanned_partitions":summary["scanned_partitions"]})
        return
    controls=json.loads((OUT/"controls.json").read_text())
    q,d,names=model_training(controls)
    diagnostics_result=diagnostics(q,d)
    save(OUT/"probability_model.json",diagnostics_result)
    valid=diagnostics_result["model_status"]=="UNVALIDATED_CONDITIONAL_TWO_CLASS_MODEL"
    source=OUT/"ranked_survivors.jsonl.gz"
    rows=[]
    with gzip.open(source,"rt",encoding="utf-8") as f:
        for line in f:
            rows.append(json.loads(line))
    for row in rows:
        y=good(row.get("Y_minus_J_AB"))
        jh=good(row.get("J_minus_H_AB"))
        if valid and y is not None and jh is not None:
            distance=min_training_distance((y,jh),q,d)
            row["nearest_labeled_colour_distance_kde_bandwidths"]=round(distance,3)
            if distance > 4:
                row["p_quasar_vs_dwarf"]=None
                row["p_quasar_vs_dwarf_prior_1pct"]=None
                row["probability_status"]="OUTSIDE_EMPIRICAL_COLOUR_SUPPORT_NO_NUMERIC_PROBABILITY"
            else:
                p1=p_q_given_q_or_d((y,jh),q,d)
                p10=p_q_given_q_or_d((y,jh),q,d,.01)
                if p1 is not None and p10 is not None:
                    row["p_quasar_vs_dwarf"]=round(p1,7)
                    row["p_quasar_vs_dwarf_prior_1pct"]=round(p10,7)
                    row["probability_status"]="UNVALIDATED_CONDITIONAL_Q_VS_DWARF_PRIOR_0P001"
                else:
                    row["p_quasar_vs_dwarf"]=None
                    row["p_quasar_vs_dwarf_prior_1pct"]=None
                    row["probability_status"]="ZERO_EMPIRICAL_DENSITY_NO_NUMERIC_PROBABILITY"
        else:
            row["p_quasar_vs_dwarf"]=None
            row["p_quasar_vs_dwarf_prior_1pct"]=None
            row["probability_status"]="UNAVAILABLE_INSUFFICIENT_EMPIRICAL_CONTROLS_OR_COLOURS"
    # Rank on fixed photometric heuristic, NOT model output; a two-class
    # posterior cannot be used as a discovery-quality score.
    ranked=sorted(rows,key=lambda r:(-r["rank_score"],r["id"]))
    final=OUT/"ranked_vetted.jsonl.gz"
    with gzip.open(final,"wt",encoding="utf-8",compresslevel=6) as f:
        for row in ranked:
            f.write(json.dumps(row,sort_keys=True,allow_nan=False)+"\n")
    save(OUT/"top_candidates.json",ranked[:100])
    supported=sorted((r for r in ranked if r["p_quasar_vs_dwarf"] is not None),
                     key=lambda r:(-r["p_quasar_vs_dwarf"],-r["rank_score"],r["id"]))
    save(OUT/"probability_supported_candidates.json",supported[:100])
    external=[]
    work=list(enumerate(ranked[:TOP_EXTERNAL]))
    with ThreadPoolExecutor(max_workers=4) as pool:
        fut={pool.submit(vet_one,row,i):i for i,row in work}
        for f in as_completed(fut):
            i=fut[f]
            try:
                external.append((i,f.result()))
            except Exception as ex:
                external.append((i,{"id":work[i][1]["id"],
                           "status":"EXTERNAL_VETTING_EXCEPTION",
                           "error":type(ex).__name__+": "+str(ex)[:500]}))
            print("EXTERNAL_VET",i,"/",len(work),flush=True)
    external=[r for _,r in sorted(external)]
    save(OUT/"external_receipts.json",external)
    by_id={r["id"]:r for r in external}
    final_top=[]
    for r in ranked[:100]:
        x=dict(r)
        match=by_id.get(r["id"])
        x["external_vetting"]=match
        x["image_level_artifact_status"]=(
           match.get("queries",{}).get("Euclid_image_metadata",{}).get("image_level_artifact_status",
                 "NOT_EVALUATED") if match else "NOT_EVALUATED")
        x["followup_priority"]=(
            "REJECT_POSSIBLE_FOREGROUND" if match and
            match.get("diagnostics",{}).get("foreground_motion_diagnostic","").startswith("REJECT")
            else "UNVETTED_NO_SPECTROSCOPIC_CONFIRMATION")
        final_top.append(x)
    save(OUT/"top_candidates_with_checks.json",final_top)
    status={"status":"PHOTOMETRIC_CANDIDATE_RANKING_ONLY",
            "survey_scan_complete":True,
            "full_matched_control_training_Q":len(q),
            "full_matched_control_training_dwarf":len(d),
            "probability_model_status":diagnostics_result["model_status"],
            "total_ranked_survivors":len(ranked),
            "external_records_attempted":len(external),
            "external_records_with_any_error":sum(
                 bool(x.get("status")=="EXTERNAL_VETTING_EXCEPTION" or
                      any(y.get("status")=="QUERY_ERROR" for y in x.get("queries",{}).values()))
                 for x in external),
            "image_level_artifact_clear_count":0,
            "spectroscopically_confirmed_new_quasars":0,
            "all_spectra_and_literature_complete":False,
            "full_known_z_gt_6p5_footprint_recovery_fraction":None,
            "posterior_is_calibrated":False,
            "posterior_out_of_empirical_colour_support_count":sum(
                 x["probability_status"]=="OUTSIDE_EMPIRICAL_COLOUR_SUPPORT_NO_NUMERIC_PROBABILITY"
                 for x in ranked),
            "conditional_probabilities_available_count":sum(
                 x["p_quasar_vs_dwarf"] is not None for x in ranked),
            "discovery_level":False,
            "run_utc":utc()}
    save(OUT/"vetting_status.json",status)
    print("VETTING",json.dumps(status),flush=True)

if __name__=="__main__":
    main()
