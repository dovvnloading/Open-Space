#!/usr/bin/env python3
"""Gaia DR3 astrometric-stellar filter on the post-SDSS 4XMM host cohort.

Galaxy identity and complete blazar rejection are not implied by this check.
Failed Gaia batches remain UNKNOWN, never certified star-free.
"""
import csv
import datetime as dt
import hashlib
import io
import json
import math
import time
from pathlib import Path
from xml.sax.saxutils import escape

import requests
BASE=Path(__file__).resolve().parent
OUT=BASE/"outputs"
SOURCE=OUT/"archive_sdss_after_exclusions.csv"
FROZEN=OUT/"archive_sdss_stellar_blazar_audit.json"
S=requests.Session()
S.headers["User-Agent"]="Open-Space-XMM-QPE-Gaia-DR3-stellar-audit"
ENDPOINTS=[
 "https://gea.esac.esa.int/tap-server/tap/sync",
 "https://gaia.ari.uni-heidelberg.de/tap/sync"
]
ADQL=("SELECT x.detid, g.source_id, g.ra AS gaia_ra, g.dec AS gaia_dec,"
      " g.parallax, g.parallax_over_error, g.pmra, g.pmra_error,"
      " g.pmdec, g.pmdec_error, g.ruwe "
      " FROM TAP_UPLOAD.xray AS x JOIN gaiadr3.gaia_source AS g"
      " ON 1=CONTAINS(POINT('ICRS',g.ra,g.dec),"
      " CIRCLE('ICRS',x.ra,x.dec,0.000277777778))")
HEADER=("detid","source_id","gaia_ra","gaia_dec","parallax",
        "parallax_over_error","pmra","pmra_error","pmdec","pmdec_error","ruwe")

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024**2),b""):h.update(block)
    return h.hexdigest()

def xml_upload(rows):
    cells="".join("<TR><TD>%s</TD><TD>%.9f</TD><TD>%.9f</TD></TR>" % (
        escape(r["detid"]),float(r["ra"]),float(r["dec"])) for r in rows)
    return ("""<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3" version="1.3">
<RESOURCE type="results"><TABLE name="xray">
<FIELD name="detid" datatype="char" arraysize="*"/>
<FIELD name="ra" datatype="double" unit="deg"/>
<FIELD name="dec" datatype="double" unit="deg"/>
<DATA><TABLEDATA>""" + cells +
        """</TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>""").encode()

def get_gaia(group):
    payload=xml_upload(group)
    errors=[]
    for endpoint in ENDPOINTS:
        for attempt in (1,2):
            t0=time.monotonic()
            try:
                with S.post(endpoint,data={"REQUEST":"doQuery","LANG":"ADQL",
                       "FORMAT":"csv","QUERY":ADQL,"UPLOAD":"xray,param:xray"},
                       files={"xray":("xray.xml",payload,"application/x-votable+xml")},
                       timeout=(12,80)) as resp:
                    body=resp.content.decode("utf-8","replace")
                    resp.raise_for_status()
                    lines=[x for x in body.splitlines() if x and not x.startswith("#")]
                    if not lines or not lines[0].lower().startswith("detid,"):
                        raise ValueError("TAP non-CSV response: "+body[:450])
                    data=list(csv.DictReader(io.StringIO("\n".join(lines))))
                    fields=set(lines[0].lower().split(","))
                    if not {"detid","source_id","parallax","pmra","pmdec"}.issubset(fields):
                        raise ValueError("TAP columns missing: "+lines[0])
                    return data,{"ok":True,"endpoint":endpoint,"attempt":attempt,
                      "query_sha256":hashlib.sha256(ADQL.encode()).hexdigest(),
                      "n_targets":len(group),"n_rows":len(data),
                      "bytes":len(body),"elapsed_s":round(time.monotonic()-t0,2)}
            except Exception as exc:
                errors.append({"endpoint":endpoint,"attempt":attempt,
                    "error":repr(exc)[:500],"elapsed_s":round(time.monotonic()-t0,2)})
    return None,{"ok":False,"errors":errors,"n_targets":len(group)}

def finite(v):
    try:
        z=float(v)
        return z if math.isfinite(z) else None
    except (TypeError,ValueError):return None

def stellar(r):
    z=finite(r.get("parallax_over_error"))
    if z is not None and z>=5 and (finite(r.get("parallax")) or 0)>0:
        return True,"Gaia_DR3_positive_parallax_snr_ge_5"
    for ax in ("pmra","pmdec"):
        val,err=finite(r.get(ax)),finite(r.get(ax+"_error"))
        if val is not None and err is not None and err>0 and abs(val)/err>=5:
            return True,"Gaia_DR3_"+ax+"_snr_ge_5"
    return False,"no_significant_Gaia_astrometric_motion"

def main():
    cfg=json.loads((BASE/"control_batch.json").read_text())
    assert cfg.get("stage")=="archive_gaia"
    sdss=json.loads(FROZEN.read_text())
    assert sdss["status"]=="SDSS_SUBSET_EXCLUSIONS_COMPLETE"
    rows=list(csv.DictReader(SOURCE.open(newline="")))
    assert len(rows)==sdss["exclusion_counts"]["n_remaining_provisional"]==989
    assert len({r["detid"] for r in rows})==989
    start=int(cfg["start"]);stop=int(cfg["stop"]);size=int(cfg["chunk_size"])
    assert 0<=start<stop<=989 and 1<=size<=50
    subset=rows[start:stop]
    tag=f"{start:05d}_{stop:05d}"
    output=OUT/f"archive_gaia_crossmatch_{tag}.json"
    report={"utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "status":"running","catalogue":"Gaia DR3 gaia_source, positional separation <=1arcsec",
            "criteria":"Positive parallax SNR >=5 or absolute single component proper-motion SNR >=5",
            "input_manifest_sha256":sha(SOURCE),"start":start,"stop":stop,
            "total_post_sdss":989,"target_count":len(subset),
            "queries":[],"source_rows":[],"stellar_flags":[],
            "missing_chunks":[],"archive_screen_permitted":False,
            "limitations":["Gaia stars are potential line-of-sight matches; host verification required",
               "No comprehensive SIMBAD/4LAC/blazar matching",
               "No Legacy DR10 full-sky galaxy sample",
               "No new-source EPIC-pn light curve screened"]}
    for i in range(0,len(subset),size):
        chunk=subset[i:i+size]
        data,meta=get_gaia(chunk)
        meta.update(index_start=start+i,index_stop=start+i+len(chunk))
        report["queries"].append(meta)
        if data is None:
            report["missing_chunks"].append(start+i)
        else:
            ids={r["detid"] for r in chunk}
            for r in data:
                if r.get("detid") not in ids:continue
                report["source_rows"].append(r)
                flag,why=stellar(r)
                if flag:
                    report["stellar_flags"].append({"detid":r["detid"],
                         "gaia_source_id":r["source_id"],"reason":why,
                         "parallax":r.get("parallax"),"parallax_over_error":r.get("parallax_over_error"),
                         "pmra":r.get("pmra"),"pmdec":r.get("pmdec"),
                         "ruwe":r.get("ruwe")})
        output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    report["status"]="COMPLETE" if not report["missing_chunks"] else "INCOMPLETE"
    report["n_gaia_matches"]=len(report["source_rows"])
    report["n_stellar_flagged_detections"]=len({r["detid"] for r in report["stellar_flags"]})
    report["n_observations_screened_pn"]=0
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:report[k] for k in ("status","target_count","n_gaia_matches",
      "n_stellar_flagged_detections","missing_chunks")},indent=2),flush=True)

if __name__=="__main__":main()
