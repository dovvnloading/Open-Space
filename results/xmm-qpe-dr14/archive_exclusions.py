#!/usr/bin/env python3
"""Conservative SDSS stellar/QSO screening and ESA PPS access diagnostics.
Only *verified* exclusions are made; missing Gaia/SIMBAD coverage is a blocker.
"""
import csv
import datetime as dt
import hashlib
import io
import json
import time
from pathlib import Path
from urllib.parse import urlencode

import requests

BASE=Path(__file__).resolve().parent
OUT=BASE/"outputs"
SOURCE=OUT/"archive_pn_point_freeze.csv"
FROZEN=OUT/"archive_pn_point_freeze.json"
SKY="https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
SESSION=requests.Session()
SESSION.headers.update({"User-Agent":"Open-Space/XMM-QPE-host-exclusions-audit"})
BATCH=70

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()

def sql_rows(query):
    for attempt in range(1,4):
        try:
            r=SESSION.get(SKY,params={"cmd":query,"format":"csv"},timeout=(12,55))
            r.raise_for_status()
            clean="\n".join(x for x in r.text.splitlines() if not x.startswith("#"))
            if not clean.startswith(("detid,", "bestObjID,", "bestobjid,")):
                raise ValueError("Unexpected SQL response "+r.text[:280])
            return list(csv.DictReader(io.StringIO(clean))), {
                 "ok":True,"attempt":attempt,"bytes":len(r.content),
                 "sha256_sql":hashlib.sha256(query.encode()).hexdigest()}
        except Exception as exc:
            problem=repr(exc)[:350]
            time.sleep(attempt)
    return None,{"ok":False,"error":problem,
                 "sha256_sql":hashlib.sha256(query.encode()).hexdigest()}

def probes():
    rec=[]
    obsids=["0872392901","0306630101","0200730401","0851180401"]
    for obs in obsids:
        params={"obsno":obs,"level":"PPS","extension":"FTZ",
                "instname":"PN","name":"SRCTSR"}
        endpoints=["https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio",
                   "https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio/"]
        for endpoint in endpoints[:1]:
            u=endpoint+"?"+urlencode(params)
            t0=time.monotonic()
            entry={"obsid":obs,"url":u}
            try:
                with SESSION.get(u,stream=True,timeout=(10,22)) as response:
                    entry.update(status=response.status_code,
                                 content_type=response.headers.get("Content-Type"),
                                 content_length=response.headers.get("Content-Length"))
                    # Limit stream to first 512 bytes, do not actually download a tarball.
                    if response.status_code==200:
                        it=response.iter_content(512)
                        b=next(it,b"")
                        entry["prefix_hex"]=b[:24].hex()
                        entry["first_block_bytes"]=len(b)
            except Exception as exc:
                entry["error"]=repr(exc)[:340]
            entry["elapsed_s"]=round(time.monotonic()-t0,2)
            rec.append(entry)
    return rec

def main():
    config=json.loads((BASE/"control_batch.json").read_text())
    assert config.get("stage")=="archive_exclusions"
    freeze=json.loads(FROZEN.read_text())
    assert freeze["status"]=="PROVISIONAL_SAMPLE_FROZEN" and sha(SOURCE)==freeze["manifest_sha256"]
    rows=list(csv.DictReader(SOURCE.open(newline="")))
    assert len(rows)==1108
    info={"run_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
          "input_manifest_sha256":sha(SOURCE),"input_n":len(rows),
          "criteria":"SDSS DR18 SpecObjAll class STAR/QSO or BL-LAC subclass; clean SDSS stellar PhotoPrimary within 1 arcsec. Missing spectra -> unclassified.",
          "queries":[],"spec_records":[],"nearby_stars":[],
          "status":"RUNNING","ppsn_probe":[],"exclusion_ids":[],
          "exclusion_counts":{},"scientific_blocks":[
              "GAIA DR3 parallax/proper-motion rejection NOT established",
              "SIMBAD/4LAC complete BL Lac and blazar classification NOT established",
              "SDSS footprint is not Legacy DR10 all-sky galaxy-nucleus sample",
              "No complete EPIC-pn light curves, independent particle flare GTI, epatplot, or negative control populations"]}
    failed=[]
    for start in range(0,len(rows),BATCH):
        chunk=rows[start:start+BATCH]
        ids=",".join(r["host_objid"] for r in chunk)
        spec=("SELECT bestObjID, class, subClass, z, zWarning, sciencePrimary "
              "FROM SpecObjAll WHERE bestObjID IN ("+ids+")")
        s,log=sql_rows(spec)
        log.update(method="spectroscopy",start=start,stop=start+len(chunk))
        info["queries"].append(log)
        if s is None:failed.append(("spectroscopy",start))
        else:
            validids={r["host_objid"] for r in chunk}
            for r in s:
                if str(r.get("bestObjID","")) in validids:
                    info["spec_records"].append(r)
        values=",".join("('%s',%.8f,%.8f)"%(
              r["detid"],float(r["ra"]),float(r["dec"])) for r in chunk)
        star=("SELECT x.detid,p.objID AS star_objid,p.ra AS star_ra,"
              "p.dec AS star_dec FROM (VALUES "+values+") AS x(detid,ra,dec)"
              " CROSS APPLY dbo.fGetNearbyObjEq(x.ra,x.dec,0.01666666667) n"
              " JOIN PhotoPrimary p ON n.objID=p.objID WHERE p.type=6 AND p.clean=1 AND p.mode=1")
        a,log=sql_rows(star)
        log.update(method="nearby_stars",start=start,stop=start+len(chunk))
        info["queries"].append(log)
        if a is None:failed.append(("nearby_stars",start))
        else:
            validids={r["detid"] for r in chunk}
            info["nearby_stars"] += [x for x in a if x.get("detid") in validids]
        if (start//BATCH+1)%8==0:
            print("SDSS EXCLUSION CHUNKS",start//BATCH+1,"of",len(rows)//BATCH+1,flush=True)
    if not failed:
        host_to_det={r["host_objid"]:r["detid"] for r in rows}
        flagged={}
        for s in info["spec_records"]:
            c=(s.get("class") or "").strip().upper()
            sub=(s.get("subClass") or "").strip().upper()
            if c in ("STAR","QSO") or "BL LAC" in sub or "BLAZAR" in sub:
                det=host_to_det.get(str(s.get("bestObjID","")))
                if det:
                    flagged.setdefault(det,set()).add("sdss_spec_"+c+"_"+sub[:30])
        for s in info["nearby_stars"]:
            flagged.setdefault(s["detid"],set()).add("sdss_star_within_1arcsec")
        info["exclusion_ids"]=[{"detid":k,"reasons":sorted(v)} for k,v in sorted(flagged.items())]
        info["exclusion_counts"]={
            "unique_sdss_spectroscopic_hosts":len({x["bestObjID"] for x in info["spec_records"]}),
            "n_sdss_spectroscopic_rows":len(info["spec_records"]),
            "n_nearby_sdss_stars":len(info["nearby_stars"]),
            "n_excluded_confirmed_stellar_or_broadline_agn":len(flagged),
            "n_remaining_provisional":len(rows)-len(flagged)}
    info["failed_chunks"]=failed
    info["status"]="SDSS_PARTIAL" if failed else "SDSS_SUBSET_EXCLUSIONS_COMPLETE"
    info["ppsn_probe"]=probes()
    OUT.joinpath("archive_sdss_stellar_blazar_audit.json").write_text(
        json.dumps(info,indent=2,allow_nan=False)+"\n")
    with (OUT/"archive_sdss_after_exclusions.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        if not failed:
            rejects={x["detid"] for x in info["exclusion_ids"]}
            w.writerows(r for r in rows if r["detid"] not in rejects)
    print(json.dumps({"status":info["status"],"exclusion_counts":info["exclusion_counts"],
        "failed":failed[:5],"ppsn_probe":[{"obsid":p["obsid"],"status":p.get("status"),
             "error":p.get("error"),"first_block_bytes":p.get("first_block_bytes")}
               for p in info["ppsn_probe"]]},indent=2),flush=True)

if __name__=="__main__":
    main()
