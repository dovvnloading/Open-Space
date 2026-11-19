#!/usr/bin/env python3
"""Audited SDSS-DR18 galaxy-centroid crossmatch; NOT a QPE detection screen.

This deliberately records an SDSS-footprint subset rather than mislabelling
all 4XMM-DR14 X-ray detections as galaxy nuclei. No hidden sampled denominator.
"""
import csv
import datetime as dt
import hashlib
import io
import json
import math
import os
from pathlib import Path
import time

import requests

BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
API = requests.Session()
API.headers.update({"User-Agent": "Open-Space/4XMM-DR14-nuclear-host-audit-1.0"})
FIELDS = ["detid", "xray_ra", "xray_dec", "obs_id", "srcid", "host_objid",
          "host_ra", "host_dec", "host_type", "host_mode", "host_clean",
          "modelmag_r", "petror50_r", "offset_arcsec"]
SQL_SIGNATURE = ("SDSS DR18 PhotoPrimary type=3 (galaxy), mode=1, clean=1; "
                 "dbo.fGetNearbyObjEq at <=2.0 arcsec, independently "
                 "recomputed great-circle separation <=2 arcsec.")

def sha(path):
    dig = hashlib.sha256()
    with path.open("rb") as fp:
        for blk in iter(lambda: fp.read(2**20), b""):
            dig.update(blk)
    return dig.hexdigest()

def screen_sql(rows):
    values = ",".join("('%s',%.8f,%.8f)" % (
        r["detid"], float(r["ra"]), float(r["dec"])) for r in rows)
    # IDs and coordinates are from the checksum-verified frozen local manifest.
    return ("SELECT x.detid, p.objID AS host_objid, p.ra AS host_ra, "
            "p.dec AS host_dec, p.type AS host_type, p.mode AS host_mode, "
            "p.clean AS host_clean, p.modelMag_r AS modelmag_r, "
            "p.petroR50_r AS petror50_r "
            "FROM (VALUES " + values + ") AS x(detid,ra,dec) "
            "CROSS APPLY dbo.fGetNearbyObjEq(x.ra,x.dec,0.033333333333) AS n "
            "JOIN PhotoPrimary AS p ON n.objID=p.objID "
            "WHERE p.type=3 AND p.mode=1 AND p.clean=1")

def sep(ra1, dec1, ra2, dec2):
    # Exact spherical distance (safe across RA meridian and high declination).
    dra = math.radians((ra1-ra2+180) % 360 - 180)
    dec1,dec2 = math.radians(dec1),math.radians(dec2)
    s = math.sin((dec1-dec2)/2)**2 + math.cos(dec1)*math.cos(dec2)*math.sin(dra/2)**2
    return 2*math.asin(min(1.,math.sqrt(max(0.,s))))*180/math.pi*3600

def query(rows, timeout_s=38):
    sql = screen_sql(rows)
    for attempt in range(3):
        try:
            resp=API.get(URL,params={"cmd":sql,"format":"csv"},timeout=(10,timeout_s))
            resp.raise_for_status()
            lines=[x for x in resp.text.splitlines() if not x.startswith("#")]
            if not lines or not lines[0].lower().startswith("detid,"):
                raise ValueError("non-CSV SDSS response: "+resp.text[:400])
            data=list(csv.DictReader(io.StringIO("\n".join(lines))))
            if not set(("detid","host_objid","host_ra","host_dec")).issubset(set(data[0]) if data else set(lines[0].split(","))):
                raise ValueError("SDSS fields missing: "+lines[0])
            return data,{"ok":True,"http":resp.status_code,"bytes":len(resp.content),
                         "returned":len(data),"sql_sha256":hashlib.sha256(sql.encode()).hexdigest(),
                         "attempt":attempt+1}
        except Exception as exc:
            error=repr(exc)[:450]
            time.sleep(min(1+attempt,3))
    return None,{"ok":False,"error":error,"sql_sha256":hashlib.sha256(sql.encode()).hexdigest()}

def main():
    config=json.loads((BASE/"control_batch.json").read_text())
    assert config.get("stage")=="archive_hosts"
    frozen=json.loads((OUT/"frozen_selection.json").read_text())
    source=OUT/"frozen_detections.csv"
    verified=sha(source)
    assert verified==frozen["full_manifest_sha256"], "FROZEN MANIFEST HASH MISMATCH"
    with source.open(newline="") as fp:
        parents=list(csv.DictReader(fp))
    assert len(parents)==frozen["eligible_count"]==33488
    ids=[r["detid"] for r in parents]
    assert len(set(ids))==len(ids), "Duplicate detection IDs in parent"
    start=int(config.get("start",0))
    stop=int(config.get("stop",len(parents)))
    size=int(config.get("chunk_size",80))
    assert 0<=start<stop<=len(parents) and 1<=size<=100
    tag=f"{start:05d}_{stop:05d}"
    outfile=OUT/f"archive_sdss_matches_{tag}.csv"
    receipt=OUT/f"archive_sdss_receipt_{tag}.json"
    checks={"as_of_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "batch":tag,"status":"running", "catalog":"SDSS DR18 PhotoPrimary",
            "parent_n":len(parents), "input_manifest_sha256":verified,
            "thresholds_sha256":sha(BASE/"control_thresholds.json"),
            "criteria":SQL_SIGNATURE,"start":start,"stop":stop,
            "target_count":stop-start,"api":URL,"chunk_size":size,
            "queries":[],"n_matched":0,"n_matches":0,
            "archive_screen_permitted":False,
            "limitations":["SDSS imaging footprint only; not a full LS DR10 or all-sky sample.",
                "The SDSS morphological galaxy tag and 2 arcsec centroids are provisional; "
                "4XMM EP_EXTENT/PN exposure, Gaia/star, known blazars and multiple matches remain to be vetoed.",
                "No EPIC-pn light curves or candidate vetting was performed in this stage."]}
    seen=set()
    missing=[]
    t0=time.monotonic()
    def store_report():
        checks["n_matched"]=len(seen)
        checks["n_queries"]=len(checks["queries"])
        checks["matched_file_sha256"]=sha(outfile)
        receipt.write_text(json.dumps(checks,indent=2,allow_nan=False)+"\n")
    with outfile.open("w",newline="") as fp:
        wr=csv.DictWriter(fp,fieldnames=FIELDS);wr.writeheader()
        for first in range(start,stop,size):
            if time.monotonic()-t0>42*60:
                missing.extend(list(range(first,stop,size)))
                break
            batch=parents[first:min(first+size,stop)]
            data,metadata=query(batch)
            metadata.update(first=first, last=first+len(batch)-1)
            checks["queries"].append(metadata)
            if data is None:
                missing.append(first)
                continue
            lookup={r["detid"]:r for r in batch}
            for d in data:
                x=lookup.get(str(d.get("detid","")).strip())
                if x is None:
                    missing.append(first)
                    metadata["unexpected_id"]=d.get("detid")
                    continue
                try:
                    host_ra=float(d["host_ra"]);host_dec=float(d["host_dec"])
                    offset=sep(float(x["ra"]),float(x["dec"]),host_ra,host_dec)
                    if not math.isfinite(offset) or offset>2.001:
                        metadata["invalid_sep"]=metadata.get("invalid_sep",0)+1
                        continue
                    row={"detid":x["detid"],"xray_ra":x["ra"],"xray_dec":x["dec"],
                         "obs_id":x["obs_id"],"srcid":x["srcid"],
                         "host_objid":str(d["host_objid"]),"host_ra":host_ra,
                         "host_dec":host_dec,"host_type":d["host_type"],
                         "host_mode":d["host_mode"],"host_clean":d["host_clean"],
                         "modelmag_r":d.get("modelmag_r"),
                         "petror50_r":d.get("petror50_r"),
                         "offset_arcsec":round(offset,4)}
                    wr.writerow(row)
                    seen.add(x["detid"])
                    checks["n_matches"]+=1
                except Exception as exc:
                    metadata["parse_error"]=repr(exc)[:300]
                    missing.append(first)
            fp.flush()
            if len(checks["queries"])%20==0:
                store_report()
                print("SDSS HOST CHUNKS",len(checks["queries"]),
                      "of",math.ceil((stop-start)/size),"matches",len(seen),flush=True)
        checks["failed_chunk_starts"]=sorted(set(missing))
        checks["queries_completed"]=sum(int(q["ok"]) for q in checks["queries"])
        checks["status"]="batch_complete" if not missing else "PARTIAL_NOT_FROZEN"
        checks["completed_input_count"]=sum(
            q["last"]-q["first"]+1 for q in checks["queries"] if q["ok"])
        store_report()
    print(json.dumps({k:checks[k] for k in (
        "status","target_count","completed_input_count","n_matched","n_matches",
        "failed_chunk_starts","matched_file_sha256")},indent=2),flush=True)

if __name__=="__main__":
    main()
