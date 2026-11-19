#!/usr/bin/env python3
"""Trigger-only crosschecks against frozen Gaia/SDSS and public 4LAC/SIMBAD.

Nonmatches are never evidence of no blazar or no transient unless the query
returned a valid, complete response and its search radius is stated. A LAT
gamma-ray position is NOT an optical counterpart association.
"""
import csv, datetime as dt, hashlib, io, json, time
from pathlib import Path
import numpy as np
import requests
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u

from heasarc_batch import OUT, SESSION, manifest, save

FERMI="https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4LACDR3/"
SIMBAD="https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
MAX_SIZE=20_000_000
def utc():return dt.datetime.now(dt.timezone.utc).isoformat()

def fetch(url,params=None,attempts=3):
    rec={"url":url,"params":params,"attempts":[]}
    for n in range(1,attempts+1):
        t=time.monotonic()
        try:
            with SESSION.get(url,params=params,timeout=(12,30),stream=True) as resp:
                status=resp.status_code
                resp.raise_for_status()
                if int(resp.headers.get("Content-Length","0") or "0")>MAX_SIZE:
                    raise ValueError("20MB catalog response ceiling")
                chunks=[];length=0; h=hashlib.sha256()
                for chunk in resp.iter_content(2**16):
                    if chunk:
                        length+=len(chunk)
                        if length>MAX_SIZE:raise ValueError("stream exceeds 20MB ceiling")
                        chunks.append(chunk);h.update(chunk)
                data=b"".join(chunks)
            rec["attempts"].append({"attempt":n,"http_status":status,"bytes":length,
                                  "sha256":h.hexdigest(),"elapsed_s":round(time.monotonic()-t,2),
                                  "status":"verified_http_transfer"})
            rec["status"]="success"
            return data,rec
        except Exception as exc:
            rec["attempts"].append({"attempt":n,"elapsed_s":round(time.monotonic()-t,2),
                                   "error":repr(exc)[:240]})
            if n<attempts:time.sleep(n)
    rec["status"]="failed"
    return None,rec

def angular_arcsec(ra,dec,catra,catdec):
    return SkyCoord(float(ra)*u.deg,float(dec)*u.deg,frame="icrs").separation(
        SkyCoord(float(catra)*u.deg,float(catdec)*u.deg,frame="icrs")).arcsec

def strval(row,key):
    if key not in row.array.names:return None
    a=row[key]
    if isinstance(a,bytes):return a.decode("utf-8","replace").strip()
    return str(a).strip()

def lac_tables(targets):
    cat=[]; receipts=[]
    for suf in ("h","l"):
        url=FERMI+f"table-4LAC-DR3-{suf}.fits"
        raw,rec=fetch(url)
        receipts.append(rec)
        if raw is None:continue
        try:
            with fits.open(io.BytesIO(raw),memmap=False) as hd:
                data=next(x.data for x in hd if x.data is not None and
                         getattr(x.data,"names",None) and
                         {"RAJ2000","DEJ2000"}.issubset(set(x.data.names)))
                cols=list(data.names)
                for item in data:
                    ra=float(item["RAJ2000"]);dec=float(item["DEJ2000"])
                    d=min(angular_arcsec(ra,dec,x["ra"],x["dec"]) for x in targets)
                    if d>1800:continue
                    row={k:strval(item,k) for k in (
                        "Source_Name","CLASS","CLASS1","ASSOC1","ASSOC_FGL",
                        "RA_Counterpart","DEC_Counterpart","Counterpart",
                        "RAJ2000","DEJ2000","Conf_95_SemiMajor","Conf_95_SemiMinor")}
                    row["high_or_low"]=suf
                    row["gamma_position_deg"]=[ra,dec]
                    for rak,dek in (("RA_Counterpart","DEC_Counterpart"),("RA_CTP","DEC_CTP"),
                                    ("RA_Opt","DEC_Opt")):
                        if rak in cols and dek in cols:
                            try:
                                x=float(item[rak]);y=float(item[dek])
                                if np.isfinite(x) and np.isfinite(y):
                                    row["counterpart_position_deg"]=[x,y]
                            except (TypeError,ValueError):pass
                    cat.append(row)
            rec["parse_status"]="parsed_FITS_as_4LAC_DR3"
            rec["matched_rows_within_30arcmin_of_any_trigger"]=len(cat)
        except Exception as exc:
            rec["parse_status"]="FITS_PARSE_FAILED"
            rec["parse_error"]=repr(exc)[:300]
    return cat,receipts

def simbad_cone(ra,dec):
    radius_deg=5/3600
    query=(f"SELECT TOP 50 main_id, ra, dec, otype FROM basic WHERE "
           f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},{radius_deg:.9f}))=1")
    body,rec=fetch(SIMBAD,{"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"json","QUERY":query},2)
    rec["adql"]=query
    if body is None:return {"status":"query_failed","receipt":rec}
    try:
        response=json.loads(body)
        if isinstance(response,dict) and "data" in response:
            info=response.get("metadata",[])
            names=[m.get("name") for m in info]
            records=[]
            for line in response.get("data",[])[:50]:
                row=dict(zip(names,line))
                match={"name":row.get("main_id"),"otype":row.get("otype"),
                       "ra_deg":row.get("ra"),"dec_deg":row.get("dec")}
                if isinstance(row.get("ra"),(int,float)):
                    match["offset_arcsec"]=round(float(angular_arcsec(ra,dec,row["ra"],row["dec"])),3)
                records.append(match)
            return {"status":"complete_radius_5arcsec",
                    "matches":records,"receipt":rec}
        return {"status":"invalid_json_schema","receipt":rec}
    except Exception as exc:
        return {"status":"invalid_response","error":repr(exc)[:200],"receipt":rec}

def main():
    path=OUT/"heasarc_postrescue_referee.json"
    if not path.exists():path=OUT/"heasarc_aggregate_referee.json"
    summary=json.loads(path.read_text())
    flags=summary.get("triggers",summary.get("trigger_details",[]))
    cohort,fr=manifest()
    lookup={r["detid"]:r for rows in cohort.values() for r in rows}
    gaia=json.loads((OUT/"archive_gaia_crossmatch_00000_00989.json").read_text())
    sdss=json.loads((OUT/"archive_sdss_stellar_blazar_audit.json").read_text())
    spec={str(x["bestObjID"]):x for x in sdss.get("spec_records",[])}
    targets=[]
    for f in flags:
        oid=str(f["obs_id"]);did=str(f["detid"]);r=lookup[did]
        if r["obs_id"]!=oid:raise ValueError("trigger not in frozen host manifest")
        targets.append({"obs_id":oid,"detid":did,"ra":float(r["ra"]),"dec":float(r["dec"]),
                        "host_objid":r["host_objid"],"host_spectrum":spec.get(r["host_objid"]),
                        "gaia_matches":[x for x in gaia.get("source_rows",[]) if x["detid"]==did],
                        "gaia_significant_star_flag":any(x["detid"]==did for x in gaia.get("stellar_flags",[]))})
    cat,receipts=lac_tables(targets)
    for t in targets:
        t["simbad"]=simbad_cone(t["ra"],t["dec"])
        near=[]
        for entry in cat:
            gl=entry["gamma_position_deg"]
            deg=angular_arcsec(t["ra"],t["dec"],gl[0],gl[1])/3600.
            if deg>0.15:continue
            match={k:v for k,v in entry.items() if k not in ("gamma_position_deg",)}
            match["gamma_sep_deg"]=round(float(deg),5)
            if "counterpart_position_deg" in entry:
                cp=entry["counterpart_position_deg"]
                sep=angular_arcsec(t["ra"],t["dec"],cp[0],cp[1])
                match["counterpart_sep_arcsec"]=round(float(sep),2)
                match["positionally_matches_counterpart_within_5arcsec"]=bool(sep<=5)
            else:
                match["positionally_matches_counterpart_within_5arcsec"]=False
                match["association_warning"]="No precise counterpart RA/Dec; gamma position alone not proof"
            near.append(match)
        t["fermi_4LAC_DR3_within_0p15deg"]=sorted(near,key=lambda x:x["gamma_sep_deg"])
        t["fermi_query_complete"]=all(x.get("parse_status")=="parsed_FITS_as_4LAC_DR3" for x in receipts)
        matches=[m for m in near if m.get("positionally_matches_counterpart_within_5arcsec") is True]
        t["fermi_counterpart_associations_within_5arcsec"]=matches
        t["suspected_4LAC_blazar_or_AGN_association"]=bool(matches)
        # 4LAC includes radio galaxies and nonblazar AGN. A coordinate match
        # is an association requiring identification review, not a definitive
        # blazar label; no match does not rule out a radio-loud AGN.
        t["blazar_exclusion_confirmed"]=False
        if t["detid"]=="101445701010002":
            t["prior_known_source"]="XMMU J134736.6+173403, published radio-loud Seyfert-2 AGN; published 2003 XMM, 2008 Chandra and Swift QPO; particle flare at second apparent QPE peak"
        elif abs(t["ra"]-195.5006)<.004 and abs(t["dec"]-27.7827)<.004:
            t["prior_known_source"]="RX J1301.9+2747, published known QPE source (several XMM and Chandra epochs), NOT NEW"
        else:
            t["prior_known_source"]="Not identified by existing explicit known-source literature; novelty NOT established"
        t["status"]="crossmatches_attempted__all_science_gates_not_met"
    doc={"utc":utc(),"frozen_manifest_sha256":fr["final_csv_sha256"],
         "trigger_count":len(targets),"sources":targets,
         "fermi_catalog_receipts":receipts,
         "simbad_cone_radius_arcsec":5,"fermi_gamma_separation_radius_deg":0.15,
         "fermi_optical_counterpart_match_arcsec":5,
         "missing_independent_checks":["SAS epatplot event patterns","independent pn soft background energy-matched off-source event regions","verified flare good-time masks and MOS","event-list spatial Poisson injections and trial-adjusted false positives","full eRASS1, Swift, Chandra multi-epoch flux and spectra for unknown triggers","manual TNS and ATel novelty checks","complete Legacy footprint/blazar host census"],
         "candidate_qpe_confirmations":0,"remarkable_allowed":False,"ready_for_merge":False}
    save("heasarc_trigger_catalog_vetting.json",doc)
    print("HEASARC_CATALOG_VET",len(targets),"4LAC tables",
          sum(x.get("parse_status")=="parsed_FITS_as_4LAC_DR3" for x in receipts),
          "SIMBAD complete",sum(x["simbad"]["status"]=="complete_radius_5arcsec" for x in targets),flush=True)

if __name__=="__main__":
    main()
