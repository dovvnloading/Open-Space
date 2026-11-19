#!/usr/bin/env python3
"""Focused nine-position DJA/JADES/MAST evidence inventory. NO discovery or posterior claim.

Reads the frozen 27-group referee ledger; only the nine still-unresolved groups
are downloaded, and every independent-data limitation is recorded explicitly.
"""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import sys
import traceback

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import download_targets as down

LEDGER = HERE / "literature_referee_2026-11-19.json"
DR4_URL = "https://jades.herts.ac.uk/DR4/Combined_DR4_external_v1.2.1.fits"
DR4_README = "https://jades.herts.ac.uk/DR4/Readme_DR4_catalogues.md"
GOODS_PID = {1180,1181,1286,3215}

def clean(x):
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, bytes):
        return x.decode("utf-8", "replace").strip()
    if isinstance(x, (float, int, str, bool)) or x is None:
        return x
    return str(x)

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def json_safe(x):
    if isinstance(x, dict):
        return {str(k):json_safe(v) for k,v in x.items()}
    if isinstance(x, (list,tuple)):
        return [json_safe(v) for v in x]
    if isinstance(x,np.generic):
        return json_safe(x.item())
    if isinstance(x,float) and not np.isfinite(x):
        return None
    return x

def write(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    t = path.with_suffix(path.suffix + ".partial")
    t.write_text(json.dumps(json_safe(payload), indent=2, allow_nan=False, sort_keys=True) + "\n")
    t.replace(path)

def hdus(path):
    summary = {"sha256": down.sha256(path), "bytes": path.stat().st_size, "extensions": []}
    with fits.open(path, memmap=True) as h:
        h.verify("exception")
        for unit in h:
            shape = getattr(unit.data, "shape", None) if unit.data is not None else None
            head = unit.header
            row = {"name": unit.name, "shape": list(shape) if shape else None,
                   "extver": head.get("EXTVER"), "bunit": clean(head.get("BUNIT")),
                   "detector": clean(head.get("DETECTOR")),
                   "program": clean(head.get("PROGRAM")),
                   "nints": clean(head.get("NINTS")),
                   "exposure": clean(head.get("EXPOSURE"))}
            if unit.name == "DQ" and unit.data is not None:
                q = np.asarray(unit.data)
                row["nonzero_fraction"] = float(np.mean(q != 0))
                row["bitwise_or"] = int(np.bitwise_or.reduce(q.astype(np.uint64).ravel()))
            if unit.name == "SCI" and unit.data is not None:
                q = np.asarray(unit.data)
                row["finite_fraction"] = float(np.mean(np.isfinite(q)))
                row["median_absolute"] = float(np.nanmedian(np.abs(q)))
            if unit.name == "SLITS" and unit.data is not None:
                row["columns"] = list(unit.columns.names)
                row["slit_table_rows"] = int(len(unit.data))
                for key in ("SOURCE_ID","SLITID","MSAID","SHUTTER_ID"):
                    if key in unit.columns.names:
                        row["unique_"+key.lower()] = [clean(v) for v in np.unique(unit.data[key])[:15]]
            summary["extensions"].append(row)
        summary["has_coadd_SCI"] = "SCI" in h
        summary["has_2d_DQ"] = "DQ" in h and h["DQ"].data is not None
        summary["has_integration_level_DQ"] = False  # Individual DJA coadds cannot prove this.
    return summary

def dr4_index(path, groups):
    """Use identity + positional match, never ID alone (IDs repeat across tiers)."""
    all_matches = {str(g["group_id"]): [] for g in groups}
    with fits.open(path, memmap=True) as h:
        t = h["Obs_info"].data
        columns = list(t.names)
        for g in groups:
            pid_to_id = [(int(s.split("_")[0]), str(s.split("_")[1])) for s in g["program_source_ids"] if "_" in s and s.split("_")[0].isdigit()]
            if not pid_to_id or pid_to_id[0][0] not in GOODS_PID:
                continue
            pid, sid = pid_to_id[0]
            # The DR4 notes explicitly warn about duplicate NIRSpec_ID values.
            subset = np.where((t["PID"].astype(str) == str(pid)) &
                              (t["NIRSpec_ID"].astype(str) == sid))[0]
            ra, dec = g["coord_deg"]["ra"], g["coord_deg"]["dec"]
            for i in subset[:50]:
                r = t[i]
                try:
                    sep = SkyCoord(float(r["RA_TARG"])*u.deg,float(r["Dec_TARG"])*u.deg).separation(
                        SkyCoord(ra*u.deg,dec*u.deg)).arcsec
                except (ValueError, TypeError):
                    continue
                fields = ["Unique_ID","PID","TIER","NIRSpec_ID","NIRCam_ID_DR3",
                          "NIRCam_ID_DR5","z_Spec","z_Spec_flag","z_PRISM",
                          "z_R1000","nInt_Prism","nInt_G140M","nInt_G235M"]
                info = {k:clean(r[k]) for k in fields if k in columns}
                info["sep_arcsec"] = float(sep)
                info["identity_compatible_0p2_arcsec"] = bool(sep <= .2)
                all_matches[str(g["group_id"])].append(info)
    return all_matches

def mast_inventory(ra, dec):
    """List nearby MAST entries; no DQ validation and no downloaded raw integrations."""
    try:
        from astroquery.mast import Observations
        Observations.TIMEOUT = 35
        o = Observations.query_region(SkyCoord(ra*u.deg, dec*u.deg), radius=4*u.arcsec)
        entries = []
        for r in o[:250]:
            if str(r["obs_collection"]).upper() != "JWST":
                continue
            entries.append({k:clean(r[k]) for k in
                ["obs_id","proposal_id","instrument_name","dataproduct_type","calib_level"]
                if k in o.colnames})
            if len(entries) >= 25:
                break
        return {"status":"INVENTORY_NOT_EXPOSURE_VALIDATION", "jwst_rows_sampled":len(entries),
                "sample": entries, "query_radius_arcsec":4.0, "query_rows_total":len(o)}
    except Exception as e:
        return {"status":"MAST_INVENTORY_FAILED", "error":type(e).__name__+": "+str(e)[:500]}

def run(out, no_mast=False):
    ledger = json.loads(LEDGER.read_text())
    groups = [g for g in ledger["entries"] if g["referee_class"] == "UNRESOLVED_NOT_NEW"]
    assert len(groups) == 9 and sum(g["members"] for g in groups) == 10
    work = Path("/tmp/dja-pr197-nine-fits")
    report = {"created_utc":utc(), "source_ledger_sha256":down.sha256(LEDGER),
              "source_run_id":ledger["source_run_id"], "status":"INSPECTION_ONLY_NOT_VALIDATED",
              "literature_priority_status":"NOT_FULL_ADS_ARXIV_SURVEY_SEARCH",
              "dr4_reference":DR4_README, "dr4_url":DR4_URL, "targets":[]}
    for g in groups:
        target = {"group_id":g["group_id"],"coord_deg":g["coord_deg"],
                  "program_source_ids":g["program_source_ids"],"fits":[]}
        for s in g["individual_spectra"]:
            item = {"root":s["root"],"file":s["file"]}
            root,name,url = down.safe(item)
            dest = work/root/name
            obj = dict(item,url=url)
            try:
                download = down.download_one({"root":root,"file":name,"url":url},work,cap_mb=400)
                obj.update(download)
                obj["fits_hdus"] = hdus(dest)
            except Exception as e:
                obj["status"]="DOWNLOAD_OR_FITS_INSPECTION_FAILED"
                obj["error"]=type(e).__name__+": "+str(e)[:500]
            target["fits"].append(obj)
        report["targets"].append(target)
        write(out/"followup_probe.json",report)
    try:
        dr4 = work/"Combined_DR4_external_v1.2.1.fits"
        # Capped source download: published catalogue, not synthetic source IDs.
        import urllib.request
        if not dr4.exists():
            req=urllib.request.Request(DR4_URL,headers={"User-Agent":"Open-Space-PR197-referee/1.0"})
            with urllib.request.urlopen(req,timeout=70) as response, dr4.open("wb") as f:
                n=0
                while True:
                    chunk=response.read(1024*1024)
                    if not chunk:break
                    n+=len(chunk)
                    if n > 800*1024*1024:
                        raise RuntimeError("DR4 catalogue exceeded 800 MiB cap")
                    f.write(chunk)
        matches = dr4_index(dr4,groups)
        report["dr4"]={"status":"CATALOG_MATCHES_NOT_SOURCE_VALIDATION",
                       "bytes":dr4.stat().st_size,"sha256":down.sha256(dr4),
                       "matches_by_group":matches}
    except Exception as e:
        report["dr4"]={"status":"CATALOG_RETRIEVAL_OR_PARSE_FAILED",
                       "error":type(e).__name__+": "+str(e)[:650]}
    write(out/"followup_probe.json",report)
    if not no_mast:
        for i,g in enumerate(groups):
            report["targets"][i]["mast"] = mast_inventory(**{"ra":g["coord_deg"]["ra"],"dec":g["coord_deg"]["dec"]})
            write(out/"followup_probe.json",report)
    report["completed_utc"] = utc()
    report["observed_spec_fits"] = sum(
        f.get("status") == "RECEIVED_UNINDEPENDENTLY_VALIDATED_FITS_BYTES"
        for t in report["targets"] for f in t["fits"])
    report["native_integration_dq_verified"] = 0
    report["nircam_fluxes_independently_measured"] = 0
    report["empirically_calibrated_posterior_count"] = 0
    write(out/"followup_probe.json",report)
    print(json.dumps({"status":report["status"],
        "downloaded":report["observed_spec_fits"],"dr4_status":report["dr4"]["status"],
        "jades_matches":report["dr4"].get("matches_by_group",{})},default=str)[:9000])

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--no-mast",action="store_true")
    args=p.parse_args()
    run(args.out,args.no_mast)
