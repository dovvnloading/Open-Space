#!/usr/bin/env python3
"""Independent archive triage of frozen, beam-aware HI--Euclid overlaps.

Network access is limited to this GitHub Actions job. Nonmatches in source
catalogues are not image-based nondetections or stellar mass limits.
"""
import concurrent.futures as cf
import csv
import datetime as dt
import hashlib
import io
import json
import math
import threading
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from astropy.io.votable import parse_single_table
import requests
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
OUT=HERE/"output"
CFG=json.loads((HERE/"freeze.json").read_text())
USERAGENT={"User-Agent":"Open-Space-Euclid-DR1F-HI-beam-referee/1.0"}
TASKS={
    "euclid_wide":("https://irsa.ipac.caltech.edu/TAP/sync","euclid_dr1_mer_catalogue_wide_survey",
                   "object_id,ra,dec,flux_vis_3fwhm_aper,fluxerr_vis_3fwhm_aper,flux_j_3fwhm_aper,fluxerr_j_3fwhm_aper", "ra","dec"),
    "euclid_deep":("https://irsa.ipac.caltech.edu/TAP/sync","euclid_dr1_mer_catalogue_deep_survey",
                   "object_id,ra,dec,flux_vis_3fwhm_aper,fluxerr_vis_3fwhm_aper,flux_j_3fwhm_aper,fluxerr_j_3fwhm_aper", "ra","dec"),
    "legacy_dr10":("https://datalab.noirlab.edu/tap/sync","ls_dr10.tractor",
                   "ls_id,ra,dec,flux_g,flux_r,flux_z,maskbits","ra","dec"),
    "unwise":("https://irsa.ipac.caltech.edu/TAP/sync","unwise_2019",
                   "ra,dec,flux_1,dflux_1,flux_2,dflux_2,flags_unwise_1,flags_unwise_2","ra","dec"),
    "desi_dr1":("https://datalab.noirlab.edu/tap/sync","desi_dr1.zpix",
                   "targetid,mean_fiber_ra,mean_fiber_dec,z,zwarn,deltachi2,spectype,zcat_primary",
                   "mean_fiber_ra","mean_fiber_dec"),
}
def utc():return dt.datetime.now(dt.timezone.utc).isoformat()
def digest(x):return hashlib.sha256(x).hexdigest()
def num(v):
    try:
        x=float(v);return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None
def separation(ra,dec,rb,db):
    a,b,c,d=[math.radians(x) for x in (ra,dec,rb,db)]
    q=math.sin((b-d)/2)**2+math.cos(b)*math.cos(d)*math.sin((a-c)/2)**2
    return 7200*math.degrees(math.asin(min(1,math.sqrt(max(0,q)))))
def query(row,task):
    endpoint,table,cols,racol,deccol=TASKS[task]
    ra,dec=float(row.ra_deg),float(row.dec_deg)
    radius=float(row.beam_fwhm_arcsec)/2/3600
    # finite beam radius, not a 2" optical match; do not treat dense-field
    # SELECT TOP saturation as a complete counterpart search.
    cap=400 if task.startswith("euclid") else (150 if task=="legacy_dr10" else 100)
    # NOIRLab's TAP PostgreSQL parser does not implement ADQL POINT().
    # Query an enclosing RA/Dec box, then enforce the angular beam radius
    # with separation() below. Handle RA=0 wrap explicitly.
    dra=radius/max(.01,math.cos(math.radians(dec)))
    ra_low=(ra-dra)%360; ra_high=(ra+dra)%360
    dec_low=max(-90,dec-radius);dec_high=min(90,dec+radius)
    if ra_low<=ra_high:
        ra_clause=f"({racol} BETWEEN {ra_low:.8f} AND {ra_high:.8f})"
    else:
        ra_clause=f"({racol}>={ra_low:.8f} OR {racol}<={ra_high:.8f})"
    sql=(f"SELECT TOP {cap} {cols} FROM {table} WHERE "
         f"({deccol} BETWEEN {dec_low:.8f} AND {dec_high:.8f}) "
         f"AND {ra_clause}")
    rec={"task":task,"table":table,"url":endpoint,"query":sql,
         "query_sha256":digest(sql.encode()),"utc":utc()}
    try:
        r=requests.post(endpoint,data={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":sql},
                        headers=USERAGENT,timeout=(10,22))
        rec.update({"http_code":r.status_code,"bytes":len(r.content),
                    "response_sha256":digest(r.content)})
        r.raise_for_status()
        raw=r.content.decode("utf-8-sig","replace")
        if raw.lstrip().startswith("<"):
            if "<VOTABLE" not in raw[:1200].upper():
                raise ValueError("Non-VOTABLE XML/HTML: "+raw[:310].replace("\\n"," "))
            try:
                root=ET.fromstring(raw)
                errors=[el for el in root.iter() if el.tag.endswith("INFO")
                    and el.attrib.get("name")=="QUERY_STATUS"
                    and el.attrib.get("value")=="ERROR"]
                if errors:
                    raise ValueError("TAP query rejected: "+
                        ";".join((e.text or "").strip() for e in errors)[:400])
                vt=parse_single_table(io.BytesIO(r.content)).to_table()
                rows=[{str(k).lower():("" if np.ma.is_masked(rec[k]) else
                     (rec[k].decode("utf-8","replace") if isinstance(rec[k],bytes)
                      else str(rec[k]))) for k in vt.colnames} for rec in vt]
            except Exception as exc:
                raise ValueError("VOTABLE_PARSE_OR_QUERY_ERROR: "+str(exc)[:180]+
                                 "; XML prefix: "+raw[:310].replace("\\n"," ")) from exc
            rec["response_format"]="VOTABLE"
        else:
            rows=[{str(k).lower().strip():v for k,v in ent.items()}
                  for ent in csv.DictReader(io.StringIO(raw))]
            rec["response_format"]="CSV"
        if rows and (racol not in rows[0] or deccol not in rows[0]):
            raise ValueError("Expected coordinate columns absent")
        if not rows and not raw.strip():raise ValueError("Empty TAP response")
        rec.update({"status":"OK","rows":len(rows),"cap":cap,"capped":len(rows)>=cap,
                    "header":raw.splitlines()[0][:550] if raw else ""})
        # keep a compact inspectable representative set of nearest sources
        matched=[]
        for ent in rows:
            rra=num(ent.get(racol));rdec=num(ent.get(deccol))
            if rra is None or rdec is None:continue
            s=separation(ra,dec,rra,rdec)
            if s<=radius*3600+0.2:
                matched.append((s,ent))
        matched.sort(key=lambda x:x[0])
        rec["nearest_arcsec"]=round(matched[0][0],4) if matched else None
        nearest=[{"sep_arcsec":round(s,4),"row":{k:v for k,v in ent.items()}} for s,ent in matched[:4]]
        return rec,nearest
    except Exception as ex:
        rec.update({"status":"QUERY_ERROR","error":str(ex)[:550]})
        return rec,[]
def evaluate(row,results):
    statuses={k:results[k][0]["status"] for k in TASKS}
    spec=[]
    rv=num(row.vopt_kms_or_derived)
    for entry in results["desi_dr1"][1]:
        e=entry["row"];z=num(e.get("z"));warn=num(e.get("zwarn"));quality=num(e.get("deltachi2"))
        if z is None or rv is None:continue
        dv=299792.458*math.log((1+z)/(1+rv/299792.458)) if 1+rv/299792.458>0 and 1+z>0 else float("nan")
        entry["velocity_offset_kms"]=round(dv,2) if math.isfinite(dv) else None
        if abs(dv)<=350 and warn==0 and (quality is None or quality>=25):
            spec.append(entry)
    nearest_vis=[]
    for k in ("euclid_wide","euclid_deep"):
        near=results[k][1]
        if near:
            nearest_vis.append((k,near[0]))
    # A source with a DESI spectroscopic counterpart is NOT a dark galaxy.
    if spec:
        status="SPECTROSCOPIC_GALAXY_REDSHIFT_ASSOCIATION"
    elif any(x=="QUERY_ERROR" for x in statuses.values()):
        status="INCOMPLETE_ARCHIVE_QUERIES"
    elif nearest_vis and min(x[1]["sep_arcsec"] for x in nearest_vis)<=5:
        status="CENTRAL_EUCLID_CATALOG_OBJECT_UNCLASSIFIED"
    else:
        status="NO_CONFIRMED_REDSHIFT_COUNTERPART_NOT_STARLESS"
    return {
      "survey":row.survey,"hi_id":row.id,"ra_deg":row.ra_deg,"dec_deg":row.dec_deg,
      "beam_fwhm_arcsec":row.beam_fwhm_arcsec,"velocity_kms":rv,
      "w50_kms":num(row.w50_kms),"logMHI":num(row.logMHI),
      "footprint":row.coverage_label,"catalog_status":statuses,
      "spectroscopic_matches":spec,"disposition":status,
      "nearby_euclid_objects":nearest_vis,
      "catalog_nearest":{k:{"sep_arcsec":results[k][0].get("nearest_arcsec"),
                             "n_rows":results[k][0].get("rows"),"capped":results[k][0].get("capped")}
                             for k in TASKS},
      "pixel_VIS_NISP_sensitivity":None,"mstar_upper_limit_msun":None,
      "MHI_over_Mstar_lower_limit":None,
      "rotation_interpretation":"UNTESTED: W50 includes turbulent/confusion/inclination effects",
      "environment_interpretation":"UNTESTED: nearby HI peers not a tidal feature test",
      "RFI_and_spectral_cube_validation":"UNTESTED",
      "starless_confirmed":False
    }
def one(row):
    # Each target is a single serialized 5-catalog archive unit; target units
    # are parallelized, with no global hand-chosen counterparts.
    result={k:query(row,k) for k in TASKS}
    study=evaluate(row,result)
    receipts=[{"hi_id":row.id,**result[k][0]} for k in TASKS]
    study["nearest_catalog_entries"]={k:result[k][1][:2] for k in TASKS}
    return study,receipts
def main():
    OUT.mkdir(exist_ok=True)
    p=OUT/"footprint_hits.csv"
    if not p.exists():
        raise RuntimeError("All-row Euclid footprint table absent; cannot run optical followup")
    hits=pd.read_csv(p)
    # Explicit stratification: 16 high-mass/wide, 16 low-mass/low-W50 from
    # each survey where available. This is followup sampling, NOT all-row
    # optical validation. Strictly distinguish from full footprint census.
    targets=[]
    cap=int(CFG["followup_cap"])
    surveys=list(hits["survey"].unique())
    for s in surveys:
        sub=hits[(hits.survey==s)&(hits.euclid_coverage_code==2)]
        if sub.empty:continue
        budget=cap//max(1,len(surveys))
        high=sub.sort_values("logMHI",ascending=False,na_position="last").head(budget//2)
        faint=sub.sort_values(["w50_kms","logMHI"],ascending=True,na_position="last").head(budget-budget//2)
        targets.append(pd.concat((high,faint)).drop_duplicates(["survey","id"]))
    if targets:
        targ=pd.concat(targets).drop_duplicates(["survey","id"]).head(cap)
    else:targ=hits.head(0)
    print("Frozen optical-followup sample:",len(targ),[(x,int((targ.survey==x).sum())) for x in surveys],flush=True)
    result=[];all_receipts=[]
    with cf.ThreadPoolExecutor(max_workers=5) as pool:
        futures={pool.submit(one,row):row for row in targ.itertuples(index=False)}
        for f in cf.as_completed(futures):
            item=futures[f]
            try:
                study,recs=f.result();result.append(study);all_receipts.extend(recs)
                print("Checked",study["survey"],study["hi_id"],study["disposition"],flush=True)
            except Exception as exc:
                result.append({"survey":item.survey,"hi_id":item.id,"disposition":"PROCESSING_ERROR",
                    "error":str(exc)[:200],"starless_confirmed":False})
    result.sort(key=lambda x:(x["survey"],x["hi_id"]))
    all_receipts.sort(key=lambda x:(x["hi_id"],x["task"]))
    resultobj={"utc":utc(),"selection":"stratified diagnostic subset of MOC center hits, high logMHI and lowest W50; not all rows",
               "footprint_hits":len(hits),"targets_attempted":len(result),"results":result,
               "conclusion":"No dark galaxy classification or stellar mass upper limit from catalog nonmatches."}
    (OUT/"followup.json").write_text(json.dumps(resultobj,indent=2,allow_nan=False)+"\n")
    (OUT/"followup_receipts.json").write_text(json.dumps(all_receipts,indent=2,allow_nan=False)+"\n")
    summary={"review_utc":utc(),"all_footprint_hits":len(hits),
       "followup_attempted":len(result),
       "valid_desired_redshift_counterparts":sum(x.get("disposition")=="SPECTROSCOPIC_GALAXY_REDSHIFT_ASSOCIATION" for x in result),
       "sources_with_any_query_error":sum(x.get("disposition")=="INCOMPLETE_ARCHIVE_QUERIES" for x in result),
       "stellar_mass_upper_limits_measured":0,"confirmed_dark_galaxies":0,
       "all_rows_optical_coverage_proven":False,
       "gates_failed":["No validated, completeness-calibrated VIS/NISP surface-brightness limit",
                       "No forced SED flux or Mstar upper limits","No HI cube/RFI inspection",
                       "No environmental/tidal discrimination","No exhaustive optical association over all in-footprint HI rows"],
       "decision":"INCOMPLETE SCREEN; no discovery or null incidence inference"}
    (OUT/"referee.json").write_text(json.dumps(summary,indent=2)+"\n")
if __name__=="__main__":main()
