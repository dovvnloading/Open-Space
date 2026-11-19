#!/usr/bin/env python3
"""Bounded archive checks. Empty queries are NOT flux limits or image vetting."""
import csv
import hashlib
import io
import json
import math
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path

import requests

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
N=40
EPOCH_EUCLID=2024.25
CATWISE_EPOCH=2015.405
IRSA="https://irsa.ipac.caltech.edu/TAP/sync"
HEADERS={"User-Agent":"Open-Space-DR1F-independent-cold-TY/1.0"}

def utc():
    return datetime.now(timezone.utc).isoformat()

def save(file,obj):
    (OUT/file).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")

def num(x):
    try:
        f=float(x)
        return f if math.isfinite(f) else None
    except (ValueError,TypeError):return None

def sep(a,d,b,e):
    t=math.sin(math.radians((d-e)/2))**2
    t+=math.cos(math.radians(d))*math.cos(math.radians(e))*math.sin(math.radians((a-b)/2))**2
    return math.degrees(2*math.asin(math.sqrt(max(0,min(1,t)))))*3600

def tap(label,endpoint,sql,ra,dec,max_sep=3,fields=("ra","dec"),timeout=16):
    rec={"source":label,"endpoint":endpoint,"adql":sql,"utc":utc(),
         "status":"QUERY_ERROR","nearest_arcsec":None,"catalog_match":None}
    try:
        response=requests.get(endpoint,params={"REQUEST":"doQuery","LANG":"ADQL",
           "FORMAT":"csv","QUERY":sql},headers=HEADERS,timeout=(8,timeout))
        rec.update({"http_status":response.status_code,
                    "response_sha256":hashlib.sha256(response.content).hexdigest(),
                    "response_bytes":len(response.content)})
        response.raise_for_status()
        if response.text.lstrip().startswith("<"):
            raise ValueError("TAP returned XML/service error: "+response.text[:160])
        entries=list(csv.DictReader(io.StringIO(response.text.lstrip("\ufeff"))))
        rec["returned_rows"]=len(entries)
        rec["status"]="CATALOG_ROWS" if entries else "NO_CATALOG_ROWS_NOT_AN_IMAGE_LIMIT"
        hits=[]
        for r in entries:
            if label=="CatWISE2020":
                a,b=num(r.get("ra_pm")),num(r.get("dec_pm"))
                mx,my=num(r.get("pmra")),num(r.get("pmdec"))
                if None in (a,b,mx,my):continue
                dt=EPOCH_EUCLID-CATWISE_EPOCH
                a=(a+dt*mx/(3600*max(.05,math.cos(math.radians(b)))))%360
                b=b+dt*my/3600
            else:
                a,b=num(r.get(fields[0])),num(r.get(fields[1]))
            if a is None or b is None:continue
            distance=sep(ra,dec,a,b)
            radius=2.5 if label=="CatWISE2020" else max_sep
            if distance<=radius:hits.append((distance,r))
        if hits:
            dis,row=min(hits,key=lambda x:x[0])
            rec["nearest_arcsec"]=round(dis,4)
            rec["catalog_match"]=row
            if label=="CatWISE2020":
                rec["motion_propagated_match"]=True
    except Exception as e:
        rec["error"]=type(e).__name__+": "+str(e)[:250]
    return rec

def catwise(r):
    ra,dec=float(r["ra"]),float(r["dec"])
    # Position at 2024 may differ by 9 yr of motion from CatWISE epoch.
    # A 15 arcsec search misses very high-PM objects; record that limitation.
    sql=("SELECT TOP 100 source_name,ra,dec,ra_pm,dec_pm,pmra,pmdec,"
         "sigpmra,sigpmdec,w1mpro_pm,w2mpro_pm,w1snr_pm,w2snr_pm "
         "FROM catwise_2020 WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
         f"CIRCLE('ICRS',{ra:.8f},{dec:.8f},{15/3600:.9f}))")
    item=tap("CatWISE2020",IRSA,sql,ra,dec,max_sep=16,timeout=22)
    item["motion_propagated_match"]=bool(item.get("catalog_match"))
    item["position_only_without_motion_fit"]=False
    if item.get("catalog_match"):
        row=item["catalog_match"]
        a=num(row.get("ra_pm"));d=num(row.get("dec_pm"))
        x=num(row.get("pmra"));y=num(row.get("pmdec"))
        if None not in (a,d,x,y):
            ep=EPOCH_EUCLID-CATWISE_EPOCH
            predicted_dec=d+ep*y/3600
            predicted_ra=(a+ep*x/(3600*max(.05,math.cos(math.radians(d)))))%360
            residual=sep(ra,dec,predicted_ra,predicted_dec)
            item["predicted_2024p25_separation_arcsec"]=round(residual,4)
            item["motion_propagated_match"]=residual<=2.5
            item["nominal_catalog_position_only"]=False
            item["w1_minus_w2_vega"]=None
            w1,w2=num(row.get("w1mpro_pm")),num(row.get("w2mpro_pm"))
            sn1,sn2=num(row.get("w1snr_pm")),num(row.get("w2snr_pm"))
            if None not in (w1,w2,sn1,sn2) and sn1>=3 and sn2>=5:
                item["w1_minus_w2_vega"]=round(w1-w2,3)
            pm_sig=[]
            for v,e in (("pmra","sigpmra"),("pmdec","sigpmdec")):
                z,err=num(row.get(v)),num(row.get(e))
                if z is not None and err and err>0:pm_sig.append(abs(z/err))
            item["motion_component_max_sigma"]=round(max(pm_sig),3) if pm_sig else None
            item["w2_snr"]=sn2
    return item

def unwise(r):
    ra,dec=float(r["ra"]),float(r["dec"])
    sql=("SELECT TOP 30 ra,dec,flux_1,dflux_1,flux_2,dflux_2,"
         "fracflux_1,fracflux_2 FROM unwise_2019 "
         "WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
         f"CIRCLE('ICRS',{ra:.8f},{dec:.8f},{5/3600:.9f}))")
    item=tap("unWISE_2019",IRSA,sql,ra,dec,max_sep=3)
    item["flux_ratio_is_not_vega_colour"]=True
    return item

def legacy(r):
    ra,dec=float(r["ra"]),float(r["dec"])
    delta=.001/max(.1,math.cos(math.radians(dec)))
    sql=("SELECT TOP 30 ra,dec,flux_g,flux_r,flux_z FROM ls_dr10.tractor "
         f"WHERE ra BETWEEN {ra-delta:.8f} AND {ra+delta:.8f} "
         f"AND dec BETWEEN {dec-.001:.8f} AND {dec+.001:.8f}")
    return tap("Legacy_DR10_catalogue", "https://datalab.noirlab.edu/tap/sync",
               sql,ra,dec,max_sep=3,timeout=12)

def panstarrs(r):
    ra,dec=float(r["ra"]),float(r["dec"])
    url="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv"
    rec={"source":"PanSTARRS_DR2_mean","endpoint":url,"utc":utc(),
         "status":"QUERY_ERROR","catalog_match":None}
    try:
        response=requests.get(url,params={"ra":ra,"dec":dec,"radius":3/3600,
                            "nDetections.gte":1},timeout=(8,12),headers=HEADERS)
        rec.update({"http_status":response.status_code,
                    "response_sha256":hashlib.sha256(response.content).hexdigest(),
                    "response_bytes":len(response.content)})
        response.raise_for_status()
        rows=list(csv.DictReader(io.StringIO(response.text)))
        rec["returned_rows"]=len(rows)
        rec["status"]="CATALOG_ROWS" if rows else "NO_CATALOG_ROWS_NOT_AN_IMAGE_LIMIT"
        hits=[]
        for obj in rows:
            a,b=num(obj.get("raMean")),num(obj.get("decMean"))
            if a is not None and b is not None and sep(ra,dec,a,b)<3:
                hits.append((sep(ra,dec,a,b),obj))
        if hits:
            distance,obj=min(hits,key=lambda x:x[0])
            rec["nearest_arcsec"]=round(distance,4)
            rec["catalog_match"]={k:obj.get(k) for k in (
                "objID","raMean","decMean","gMeanPSFMag","rMeanPSFMag",
                "iMeanPSFMag","zMeanPSFMag","yMeanPSFMag",
                "gMeanPSFMagErr","zMeanPSFMagErr")}
    except Exception as e:
        rec["error"]=type(e).__name__+": "+str(e)[:250]
    return rec

def simbad(r):
    ra,dec=float(r["ra"]),float(r["dec"])
    sql=("SELECT TOP 20 main_id,ra,dec,otype FROM basic "
         "WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
         f"CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.0014))")
    attempts=[]
    for url in ("https://simbad.cfa.harvard.edu/simbad/sim-tap/sync",
                "https://simbad.idia.ac.za/simbad/sim-tap/sync"):
        item=tap("SIMBAD_mirror",url,sql,ra,dec,max_sep=5,timeout=14)
        item["mirror_literature_not_current_or_complete"]=True
        if item["status"]!="QUERY_ERROR":
            item["failed_mirrors"]=attempts
            return item
        attempts.append({"endpoint":url,"error":item.get("error")})
    return {"status":"QUERY_ERROR","source":"SIMBAD_mirrors","failed_mirrors":attempts}

def estimate(r,wise,calibration):
    out={"W1_minus_W2_vega":wise.get("w1_minus_w2_vega"),
         "Teff_estimate_K":None,"photometric_distance_pc":None,
         "Teff_and_distance_status":"NO_VALIDATED_SPECTRAL_SUBTYPE_OR_CALIBRATION"}
    match=wise.get("catalog_match") or {}
    w2=num(match.get("w2mpro_pm"))
    w2snr=num(match.get("w2snr_pm"))
    # Only *conditional* subtype envelope (T8,Y1); no colour-to-type
    # transfer is established for Euclid Y/J/H without spectra.
    anchors=calibration.get("anchors",{})
    points=[anchors[k] for k in ("T8","T9","Y0","Y1") if k in anchors]
    if w2 is not None and w2snr and w2snr>=5 and points and wise.get("motion_propagated_match"):
        ds=[]
        for entry in points:
            m=entry["median_M_W2"]
            scatter=max(.7,entry["scatter_M_W2_MAD"])
            ds.extend((10**((w2-m+5-2*scatter)/5),10**((w2-m+5+2*scatter)/5)))
        out["photometric_distance_pc"]=[round(min(ds),1),round(max(ds),1)]
        out["Teff_and_distance_status"]="CONDITIONAL_T8_TO_Y1_W2_ABSOLUTE_MAG_ENVELOPE_NOT_A_SPECTRAL_CLASSIFICATION"
        temps=[x["median_published_Teff_K"] for x in points if x.get("median_published_Teff_K")]
        if len(temps)>=2:
            out["Teff_estimate_K"]=[round(min(temps)),round(max(temps))]
    return out

def vet_one(slot,r,calibration):
    q={}
    for name,fn in (("catwise",catwise),("unwise",unwise),
                    ("legacy",legacy),("panstarrs",panstarrs),("simbad",simbad)):
        q[name]=fn(r)
    w=q["catwise"]
    motion=w.get("motion_component_max_sigma")
    color=w.get("w1_minus_w2_vega")
    motion_supported=bool(w.get("motion_propagated_match") and motion is not None and motion>=5)
    methane_supported=bool(color is not None and color>=1.8)
    known=bool(r["literature_near_epoch_matches"])
    simbad_hit=q["simbad"].get("catalog_match")
    known=known or bool(simbad_hit)
    evidence=("MOTION_PLUS_RED_WISE_PHOTOMETRIC_LEAD_NOT_CONFIRMED"
              if motion_supported and methane_supported and not known else
              "NO_JOINT_SIGNIFICANT_WISE_MOTION_AND_METHANE_EVIDENCE")
    return {"id":r["id"],"ra":r["ra"],"dec":r["dec"],"rank_slot":slot,
            "checked_utc":utc(),"known_or_simbad_association":known,
            "motion_supported":motion_supported,"WISE_methane_proxy_supported":methane_supported,
            "lead_status":evidence,
            "conditional_physical_properties":estimate(r,w,calibration),
            "query_receipts":q,"Euclid_NISP_PSF_checked":False,
            "multi_epoch_Euclid_pixel_artifacts_checked":False,
            "spectroscopic_type_confirmed":False}

def main():
    summary=json.loads((OUT/"summary.json").read_text())
    if summary["status"]!="SCAN_COMPLETE":
        raise RuntimeError("External checks only after all 768 partitions")
    shortlist=json.loads((OUT/"unvetted_shortlist.json").read_text())
    calibration=json.loads((OUT/"reference_calibration.json").read_text())
    checked=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(vet_one,i,r,calibration):i for i,r in enumerate(shortlist[:N])}
        for f in as_completed(futures):
            i=futures[f]
            try:res=f.result()
            except Exception as e:
                res={"id":shortlist[i]["id"],"rank_slot":i,
                    "lead_status":"QUERY_OR_PIPELINE_ERROR","error":repr(e)}
            checked.append(res)
            print("VET",i,len(futures),res["lead_status"],flush=True)
    checked.sort(key=lambda x:x["rank_slot"])
    save("external_receipts.json",checked)
    report=Counter()
    for r in checked:
        report["checked"]+=1
        report["known_or_simbad_associated"]+=bool(r.get("known_or_simbad_association"))
        report["possible_photometric_leads"]+=r.get("lead_status")=="MOTION_PLUS_RED_WISE_PHOTOMETRIC_LEAD_NOT_CONFIRMED"
        report["with_query_error"]+=any(x.get("status")=="QUERY_ERROR" for x in r.get("query_receipts",{}).values())
    status={"status":"NO_VALIDATED_DISCOVERY_PHOTOMETRIC_REVIEW",
        "full_hats_scan_complete":True,
        "top_n_examined_external":report["checked"],
        "new_astrometric_and_methane_photometric_leads":report["possible_photometric_leads"],
        "strongest_lead_is_spectroscopically_confirmed":False,
        "new_objects_with_validated_Teff_below_500_K":0,
        "new_objects_validated_colder_than_WISE_0855":0,
        "controls":summary["control_recovery"],
        "query_errors_count":report["with_query_error"],
        "NISP_PSF_and_pixel_artifact_checks_complete":False,
        "no_match_is_not_a_flux_upper_limit":True,
        "calibrated_completeness_available":False,
        "discovery_level":False}
    save("vetting_status.json",status)
    print("VETTING",json.dumps(status),flush=True)

if __name__=="__main__":
    main()
