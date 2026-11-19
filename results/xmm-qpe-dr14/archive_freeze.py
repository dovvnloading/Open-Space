#!/usr/bin/env python3
"""Freeze conservative POINT-LIKE, PN-qualified SDSS-DR18 host subcohort.

Not a validated complete LS/SDSS nuclear-galaxy parent; Gaia, SIMBAD/blazar
vetoes and per-observation background/pileup still require downstream review.
"""
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import shutil
import tempfile
import time

import numpy as np
from astropy.io import fits
import requests

BASE=Path(__file__).resolve().parent
OUT=BASE/"outputs"
CAT_MD5="47c6bd113cb8510db418665472c81b48"
SOURCES=["https://heasarc.gsfc.nasa.gov/FTP/xmm/data/catalogues/4XMM_DR14cat_v1.0.fits.gz",
         "https://xmmssc.irap.omp.eu/Catalogue/4XMM-DR14/4XMM_DR14cat_v1.0.fits.gz"]
API=requests.Session()
API.headers.update({"User-Agent":"Open-Space/DR14-SDSS-point-source-audit/1.0"})
PARENT=OUT/"frozen_detections.csv"
HOSTS=OUT/"archive_sdss_matches_00000_33488.csv"
NAME=OUT/"archive_pn_point_freeze.csv"
REPORT=OUT/"archive_pn_point_freeze.json"

def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as fp:
        for b in iter(lambda:fp.read(1024**2),b""):
            h.update(b)
    return h.hexdigest()

def norm_id(v):
    return v.decode("utf8","replace").strip() if isinstance(v,bytes) else str(v).strip()

def number(v,default=None):
    try:
        x=float(v)
        return x if math.isfinite(x) else default
    except (ValueError,TypeError):
        return default

def download(urls,path,records):
    for url in urls:
        r={"url":url,"start_utc":dt.datetime.now(dt.timezone.utc).isoformat()}
        start=time.monotonic()
        try:
            with API.get(url,stream=True,timeout=(25,90)) as resp:
                r["http"]=resp.status_code
                resp.raise_for_status()
                if int(resp.headers.get("content-length") or 0)>1250*2**20:
                    raise ValueError("catalogue compressed >1250MiB")
                md5=hashlib.md5()
                n=0
                with path.open("wb") as fp:
                    for block in resp.iter_content(2**20):
                        if block:
                            n+=len(block)
                            if n>1250*2**20 or time.monotonic()-start>18*60:
                                raise TimeoutError("catalogue stream budget exceeded")
                            md5.update(block)
                            fp.write(block)
                r.update(size_bytes=n,md5=md5.hexdigest())
                if md5.hexdigest()!=CAT_MD5:
                    raise ValueError("published 4XMM catalogue MD5 mismatch")
                r["verified"]=True
                records.append(r)
                return True
        except Exception as ex:
            r.update(error=repr(ex)[:400],verified=False)
            records.append(r)
            path.unlink(missing_ok=True)
    return False

def main():
    cfg=json.loads((BASE/"control_batch.json").read_text())
    assert cfg.get("stage")=="archive_freeze"
    frozen=json.loads((OUT/"frozen_selection.json").read_text())
    assert digest(PARENT)==frozen["full_manifest_sha256"]
    hosts_receipt=json.loads((OUT/"archive_sdss_receipt_00000_33488.json").read_text())
    assert hosts_receipt["status"]=="batch_complete"
    assert hosts_receipt["completed_input_count"]==33488
    assert digest(HOSTS)==hosts_receipt["matched_file_sha256"]
    with PARENT.open(newline="") as fp:
        xray={r["detid"]:r for r in csv.DictReader(fp)}
    assert len(xray)==33488
    with HOSTS.open(newline="") as fp:
        host_rows=list(csv.DictReader(fp))
    hmap=defaultdict(list)
    for row in host_rows:
        assert row["detid"] in xray
        hmap[row["detid"]].append(row)
    assert len(hmap)==hosts_receipt["n_matched"]==2510
    audit={"as_of_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
           "catalogue":"4XMM-DR14 v1.0", "catalogue_published_md5":CAT_MD5,
           "sdss":"SDSS DR18 PhotoPrimary clean=1,mode=1,type=3,host centroid <=2arcsec",
           "parent_manifest_sha256":digest(PARENT),
           "sdss_pairs_sha256":digest(HOSTS),
           "parent_xray_detections":len(xray),"sdss_matched_detections":len(hmap),
           "sdss_pairs":len(host_rows),
           "selection_pre_registered":"one clean galaxy match <=2 arcsec; EP_EXTENT <=6 arcsec and EP_EXTENT_ML<4; PN_ONTIME >=21600 s; PN_FLAG<=1; keep soft HR and TSERIES prefilter",
           "status":"running","download_receipts":[],
           "blocks":["SDSS footprint incompleteness: LS DR10 inaccessible to runner",
              "No Gaia high-PM/parallax source exclusion, no SIMBAD/4LAC known blazar clearance",
              "PN_PILEUP field is diagnostic, not an independently measured epatplot pileup veto",
              "No independent pn background region/GTI or event-level luminosity and spectra",
              "Host centers and classifications are SDSS morphology, not spectroscopically verified nuclei",
              "No EPIC-pn light curve accessed in this freeze stage"]}
    with tempfile.TemporaryDirectory() as temp:
        gz=Path(temp)/"dr14.fits.gz";out=Path(temp)/"dr14.fits"
        if not download(SOURCES,gz,audit["download_receipts"]):
            audit["status"]="CATALOGUE_UNAVAILABLE"
            REPORT.write_text(json.dumps(audit,indent=2)+"\n")
            raise SystemExit("Cannot verify full catalogue")
        with gzip.open(gz,"rb") as inp,out.open("wb") as fp:
            shutil.copyfileobj(inp,fp,1024**2)
        with fits.open(out,memmap=True) as hdul:
            det=next(x.data for x in hdul[1:] if isinstance(x,fits.BinTableHDU))
            assert len(det)==1035832,"Full DR14 row count changed"
            assert all(k in det.names for k in ("DETID","EP_EXTENT","EP_EXTENT_ML","PN_ONTIME","PN_FLAG"))
            indices={}
            for i in range(len(det)):
                identifier=norm_id(det["DETID"][i])
                if identifier in hmap:
                    indices[identifier]=i
            audit["matched_full_catalogue_rows"]=len(indices)
            assert len(indices)==len(hmap),"Missing SDSS matched 4XMM detection"
            rows=[]
            counts=Counter()
            for identifier,pairs in sorted(hmap.items()):
                d=det[indices[identifier]]
                exp=number(d["PN_ONTIME"])
                extent=number(d["EP_EXTENT"])
                extent_ml=number(d["EP_EXTENT_ML"])
                flag=number(d["PN_FLAG"])
                pileup=number(d["PN_PILEUP"]) if "PN_PILEUP" in det.names else None
                if len(pairs)!=1:
                    counts["ambiguous_galaxy"]+=1
                    continue
                if not (extent is not None and extent>=0 and extent<=6 and extent_ml is not None and extent_ml<4):
                    counts["not_pointlike"]+=1
                    continue
                if not (exp is not None and exp>=21600):
                    counts["no_long_epic_pn"]+=1
                    continue
                if not (flag is not None and flag<=1):
                    counts["pn_flag"]+=1
                    continue
                x=xray[identifier]; g=pairs[0]
                rows.append({"detid":identifier,"srcid":x["srcid"],"obs_id":x["obs_id"],
                             "pps_srcnum":x["pps_srcnum"],"ra":x["ra"],"dec":x["dec"],
                             "host_objid":g["host_objid"],"host_ra":g["host_ra"],
                             "host_dec":g["host_dec"],"offset_arcsec":g["offset_arcsec"],
                             "modelmag_r":g["modelmag_r"],
                             "pn_ontime_s":round(exp,1),"ep_extent":round(extent,2),
                             "ep_extent_ml":round(extent_ml,2),"pn_flag":int(flag),
                             "pn_pileup_diagnostic":pileup,
                             "cts":x["cts"],"hr3":x["hr3"],"hr4":x["hr4"],
                             "score":x["score"]})
            fieldnames=list(rows[0]) if rows else ["detid","srcid","obs_id"]
            with NAME.open("w",newline="") as fp:
                w=csv.DictWriter(fp,fieldnames=fieldnames)
                w.writeheader();w.writerows(rows)
            audit.update(status="PROVISIONAL_SAMPLE_FROZEN",
                         sdss_nuclear_point_pn_detections=len(rows),
                         unique_obsids=len({r["obs_id"] for r in rows}),
                         unique_xray_srcids=len({r["srcid"] for r in rows}),
                         exclusions=dict(counts),
                         manifest_sha256=digest(NAME),
                         host_list_complete_for_SDSS_subset=True,
                         archive_screen_permitted=False)
    REPORT.write_text(json.dumps(audit,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:audit[k] for k in ("status","parent_xray_detections",
        "sdss_matched_detections","sdss_nuclear_point_pn_detections",
        "unique_obsids","exclusions","manifest_sha256")},indent=2),flush=True)

if __name__=="__main__":
    main()
