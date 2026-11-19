#!/usr/bin/env python3
"""Freeze conservative POINT-LIKE, PN-qualified SDSS-DR18 host subcohort.

Not a validated complete LS/SDSS nuclear-galaxy parent; Gaia, SIMBAD/blazar
vetoes and per-observation background/pileup still require downstream review.
"""
import csv
import concurrent.futures
import os
import re
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
    """Complete checksum-pinned file via bounded concurrent HTTPS byte ranges."""
    url=urls[0]
    expected_size=1005466736  # Content-Range from dedicated 206 probe.
    part_size=16*1024*1024
    jobs=[(s,min(s+part_size,expected_size)-1)
          for s in range(0,expected_size,part_size)]
    nworkers=8
    started=time.monotonic()
    record={"url":url,"method":"HTTP Range 206; 8 parallel 16MiB chunks",
            "expected_size_bytes":expected_size,"expected_md5":CAT_MD5,
            "n_parts":len(jobs),"verified":False,"chunks":[]}
    fd=os.open(path,os.O_RDWR|os.O_CREAT|os.O_TRUNC,0o600)
    os.ftruncate(fd,expected_size)
    def get_chunk(bounds):
        first,last=bounds
        hdr={"Range":f"bytes={first}-{last}",
             "User-Agent":"Open-Space/DR14-SDSS-range-verified/1.0"}
        last_error="no attempts"
        for attempt in range(1,4):
            received=0; h=hashlib.sha256()
            try:
                if time.monotonic()-started>35*60:
                    raise TimeoutError("whole catalogue byte-range budget exceeded")
                with requests.get(url,headers=hdr,stream=True,
                                  timeout=(12,105)) as response:
                    cr=response.headers.get("Content-Range","")
                    if response.status_code!=206 or cr!=f"bytes {first}-{last}/{expected_size}":
                        raise ValueError(f"invalid HTTP Range: {response.status_code} {cr}")
                    for block in response.iter_content(1024**2):
                        if not block:continue
                        if received+len(block)>last-first+1:
                            raise ValueError("received over requested byte range")
                        written=os.pwrite(fd,block,first+received)
                        if written!=len(block):
                            raise IOError("short disk pwrite")
                        received+=len(block);h.update(block)
                    if received!=last-first+1:
                        raise IOError(f"range length {received} vs {last-first+1}")
                return {"first":first,"last":last,"sha256":h.hexdigest(),
                        "bytes":received,"attempt":attempt}
            except Exception as exc:
                last_error=repr(exc)[:250]
                time.sleep(min(attempt,3))
        raise RuntimeError(f"range {first}-{last} failed: {last_error}")
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=nworkers) as pool:
            futures=[pool.submit(get_chunk,x) for x in jobs]
            for f in concurrent.futures.as_completed(futures):
                result=f.result()
                record["chunks"].append(result)
                if len(record["chunks"])%10==0:
                    print("CATALOGUE RANGE",len(record["chunks"]),"of",len(jobs),
                          "elapsed",round(time.monotonic()-started,1),flush=True)
        os.fsync(fd)
    except Exception as exc:
        record["error"]=repr(exc)[:550]
        record["elapsed_s"]=round(time.monotonic()-started,1)
        records.append(record)
        return False
    finally:
        os.close(fd)
    m=hashlib.md5()
    with path.open("rb") as fp:
        for data in iter(lambda: fp.read(2**20),b""):
            m.update(data)
    record.update(elapsed_s=round(time.monotonic()-started,1),
                  bytes=path.stat().st_size,md5=m.hexdigest(),
                  n_received_parts=len(record["chunks"]))
    record["verified"]=bool(m.hexdigest()==CAT_MD5)
    records.append(record)
    return record["verified"]

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
           "selection_pre_registered":"one clean galaxy match <=2 arcsec; EP_EXTENT <=6 arcsec and EP_EXTENT_ML<4; PN_ONTIME >=21600 s; PN_FLAG 12-character quality mask with no serious bits 7/8/10/11/12, PN_8_CTS>=250; keep soft HR and TSERIES parent prefilter",
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
                flag=norm_id(d["PN_FLAG"]).upper()
                pncts=number(d["PN_8_CTS"])
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
                # 4XMM PN_FLAG is 12 chars F/T/- (not an integer);
                # severe quality bits 7/8/10/11/12 must be unset.
                # The parent prefilter already requires SUM_FLAG<=1.
                if not (len(flag)==12 and set(flag)<=set("FT")
                        and all(flag[j]=="F" for j in (6,7,9,10,11))):
                    counts["pn_flag_bad_or_missing"]+=1
                    continue
                if not (pncts is not None and pncts>=250):
                    counts["pn_insufficient_counts"]+=1
                    continue
                x=xray[identifier]; g=pairs[0]
                rows.append({"detid":identifier,"srcid":x["srcid"],"obs_id":x["obs_id"],
                             "pps_srcnum":x["pps_srcnum"],"ra":x["ra"],"dec":x["dec"],
                             "host_objid":g["host_objid"],"host_ra":g["host_ra"],
                             "host_dec":g["host_dec"],"offset_arcsec":g["offset_arcsec"],
                             "modelmag_r":g["modelmag_r"],
                             "pn_ontime_s":round(exp,1),"ep_extent":round(extent,2),
                             "ep_extent_ml":round(extent_ml,2),"pn_flag":flag,"pn_8_cts":round(pncts,1),
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
