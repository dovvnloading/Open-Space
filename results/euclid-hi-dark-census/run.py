#!/usr/bin/env python3
"""Frozen all-row blind H I × Euclid DR1-F footprint screen.

This is an eligibility census, not a dark-galaxy detector. All upstream
downloads/queries occur on the GitHub Actions runner. Missing data fail closed.
"""
from __future__ import annotations
import csv
import datetime as dt
import hashlib
import io
import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy import units as u
from mocpy import MOC

BASE = Path(__file__).resolve().parent
OUT = BASE / "output"
CACHE = OUT / "_raw_not_committed"
CFG = json.loads((BASE / "freeze.json").read_text())
HEADERS = {"User-Agent": "Open-Space-EuclidHI-DR1Foundation-census/1.0"}
C_KMS = 299792.458
HI_MHZ = 1420.40575177

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def finite(x):
    try:
        x = float(x)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None

def fetch_bytes(tag, url=None, endpoint=None, adql=None, expected=None):
    record = {"id": tag, "time_utc": utc(), "url": url or endpoint,
              "adql": adql, "query_sha256": sha(adql.encode()) if adql else None,
              "expected_sha256": expected}
    raw = None
    for attempt in range(3):
        try:
            if adql:
                res = requests.post(endpoint, data={"REQUEST": "doQuery", "LANG": "ADQL",
                    "FORMAT": "csv", "QUERY": adql}, headers=HEADERS, timeout=(15, 130))
            else:
                res = requests.get(url, headers=HEADERS, timeout=(15, 150))
            record["http_code"] = res.status_code
            res.raise_for_status()
            raw = res.content
            record["bytes"] = len(raw)
            record["sha256"] = sha(raw)
            record["content_type"] = res.headers.get("content-type")
            record["final_url"] = res.url
            if expected and sha(raw) != expected:
                raise ValueError("Catalog SHA256 differs from frozen reference")
            if raw.lstrip()[:6].lower().startswith((b"<html", b"<?xml")):
                raise ValueError("HTML/XML rather than expected data")
            CACHE.mkdir(parents=True, exist_ok=True)
            (CACHE / (tag + (".fits" if tag.startswith("footprint") else ".response"))).write_bytes(raw)
            record["status"] = "OK"
            return raw, record
        except Exception as exc:
            record["status"] = "ERROR"
            record["error"] = str(exc)[:250]
            if attempt < 2:
                time.sleep(2 ** attempt)
    return None, record

def parse_csv(raw, cid):
    # Cornell alpha.100 and FASHI have different case/field conventions.
    df = pd.read_csv(io.BytesIO(raw), low_memory=False)
    if len(df.columns) < 4:
        raise ValueError("Too few CSV columns")
    if cid == "ALFALFA_alpha100" and any(str(x).startswith("#") for x in df.columns):
        df = pd.read_csv(io.BytesIO(raw), low_memory=False, comment="#")
    return df

def norm(x):
    return re.sub("[^a-z0-9]", "", str(x).lower())

def column(df, names, required=False):
    d = {norm(x): x for x in df.columns}
    for n in names:
        if norm(n) in d:
            return d[norm(n)]
    if required:
        raise ValueError("Missing required column from " + str(names) + "; found " + str(list(df.columns)[:65]))
    return None

def numeric(df, names):
    c = column(df, names)
    return pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float) if c else np.full(len(df), np.nan)

def position(df, names, ra=False):
    c = column(df, names, required=True)
    vals = pd.to_numeric(df[c], errors="coerce")
    if vals.notna().mean() > .95:
        return vals.to_numpy(dtype=float), c
    # Some early source catalogues use sexagesimal fields.
    from astropy.coordinates import Angle
    unit = u.hourangle if ra else u.deg
    out = []
    for s in df[c]:
        try:
            out.append(Angle(str(s), unit=unit).degree)
        except (TypeError, ValueError):
            out.append(float("nan"))
    return np.array(out, dtype=float), c

def input_catalog(cfg):
    if cfg["type"] == "tap":
        data, receipt = fetch_bytes(cfg["id"], endpoint=cfg["endpoint"], adql=cfg["adql"])
    else:
        data, receipt = fetch_bytes(cfg["id"], url=cfg["url"], expected=cfg.get("sha256"))
    if data is None:
        return None, receipt
    try:
        df = parse_csv(data, cfg["id"])
        receipt["input_rows"] = len(df)
        receipt["columns"] = list(map(str, df.columns))
        return df, receipt
    except Exception as exc:
        receipt["status"] = "PARSE_ERROR"
        receipt["error"] = str(exc)[:550]
        return None, receipt

def normalize_catalog(df, cfg):
    cid = cfg["id"]
    if cid == "ALFALFA_alpha100":
        ras = ["HIra", "HI_RA", "RAdeg", "RAdegHI", "RAdeg_HI", "RA_HI", "RA"]
        decs = ["HIdec", "HI_DEC", "DEdeg", "DECdeg", "DecdegHI", "DECdeg_HI", "Dec_HI", "Dec"]
    else:
        ras = ["ra", "ra_deg", "ra_hi", "ra_deg_hi", "RAdeg", "ra_j2000", "pos_ra"]
        decs = ["dec", "dec_deg", "dec_hi", "dec_deg_hi", "DEdeg", "dec_j2000", "pos_dec"]
    ra, rac = position(df, ras, ra=True)
    dec, decc = position(df, decs)
    keys = ["AGC", "AGCNr", "AGCnr", "fashi_id", "name", "source_name", "id", "source_id", "name_hi"]
    key = column(df, keys)
    ids = df[key].astype(str).to_numpy() if key else np.array([str(i) for i in range(len(df))])
    # Preserve native measurements, and never guess a distance for a velocity-only source.
    v = numeric(df, ["cz", "vhelio", "Vhel", "v_sys", "vopt", "v_rad", "vel", "v_center"])
    freq = numeric(df, ["freq", "freq_mhz", "frequency"])
    fill = np.isfinite(freq) & (freq > 1000) & (freq < 1500) & ~np.isfinite(v)
    v[fill] = C_KMS * (HI_MHZ / freq[fill] - 1)
    width = numeric(df, ["W_50", "W50", "w50", "w50_kms", "w50_corrected", "W50c"])
    mass = numeric(df, ["log10Mass", "logMH", "logMHI", "log_mhi", "mass", "m_hi", "log10_mhi"])
    logmass = mass.copy()
    # catalogued mass is log10(M_sun) unless encoded as actual solar masses
    hi = np.isfinite(mass) & (mass > 1000)
    logmass[hi] = np.log10(mass[hi])
    logmass[~(np.isfinite(logmass) & (logmass > 4) & (logmass < 14))] = np.nan
    flux = numeric(df, ["S21", "HIflux", "FHI", "S_sum", "f_sum", "flux", "flux_int", "s_int", "f_sum_corr", "integrated_flux"])
    dist = numeric(df, ["distance", "Dist", "D", "D_Mpc", "dist_mpc"])
    # Only derive mass when there is a positive catalog distance and integrated flux.
    derive = ~np.isfinite(logmass) & (dist > 0) & (flux > 0)
    logmass[derive] = np.log10(2.356e5 * dist[derive]**2 * flux[derive])
    # Survey optical counterpart columns are *not* authoritative at Euclid DR1 depth.
    optical = column(df, ["OC", "optical_counterpart", "opt_id", "objname", "sdss_name", "OC_name"])
    out = pd.DataFrame({
        "survey": cid, "id": ids, "ra_deg": ra, "dec_deg": dec,
        "vopt_kms_or_derived": v, "w50_kms": width, "logMHI": logmass,
        "flux_jykms": flux, "distance_mpc": dist,
        "reported_optical": df[optical].astype(str).to_numpy() if optical else "",
        "beam_fwhm_arcsec": cfg["beam_fwhm_arcsec"],
        "source_row": np.arange(len(df), dtype=int),
    })
    meta = {"input_rows": len(df), "ra_column": rac, "dec_column": decc,
            "id_column": key, "optical_column": optical,
            "finite_positions": int((np.isfinite(ra) & np.isfinite(dec)).sum()),
            "finite_velocity": int(np.isfinite(v).sum()), "finite_W50": int(np.isfinite(width).sum()),
            "finite_logMHI": int(np.isfinite(logmass).sum())}
    return out, meta

def coverage(df, mocs):
    ra = df["ra_deg"].to_numpy(dtype=float)
    dec = df["dec_deg"].to_numpy(dtype=float)
    ok = np.isfinite(ra) & np.isfinite(dec) & (dec >= -90) & (dec <= 90)
    valid = np.flatnonzero(ok)
    status = np.zeros(len(df), dtype=np.uint8)
    radii = df["beam_fwhm_arcsec"].to_numpy(dtype=float) / 2 / 3600
    for moc in mocs:
        status[valid] = np.maximum(status[valid],
             moc.contains_lonlat(np.mod(ra[valid], 360)*u.deg, dec[valid]*u.deg).astype(np.uint8)*2)
        # Sample 16 great-circle bearings along half-power circle. Fringe contact is
        # flagged separately; it is NOT a proof of full beam coverage.
        near = valid[status[valid] == 0]
        for angle in np.arange(0, 2*np.pi, 2*np.pi/16):
            north = radii[near] * np.cos(angle)
            east = radii[near] * np.sin(angle)
            decx = np.clip(dec[near] + north, -89.9999, 89.9999)
            rax = np.mod(ra[near] + east/np.maximum(.001, np.cos(np.deg2rad(dec[near]))), 360)
            hit = moc.contains_lonlat(rax*u.deg, decx*u.deg)
            status[near[hit]] = np.maximum(status[near[hit]], 1)
    df["euclid_coverage_code"] = status
    df["coverage_label"] = np.select(
        [status == 2, status == 1], ["CENTER_IN_MOC", "BEAM_EDGE_SAMPLED"], default="OUTSIDE_OR_UNVALIDATED")
    return df

def env_counts(df):
    # In-survey nearby HI sources within 20 arcmin AND 500 km/s, descriptive only.
    from scipy.spatial import cKDTree
    x = df[df.euclid_coverage_code > 0].copy()
    if len(x) == 0:
        x["nearby_HI_20arcmin_dv500"] = []
        return x
    r = np.deg2rad(x.ra_deg.to_numpy())
    d = np.deg2rad(x.dec_deg.to_numpy())
    p = np.column_stack((np.cos(d)*np.cos(r), np.cos(d)*np.sin(r), np.sin(d)))
    tree = cKDTree(p)
    v = x.vopt_kms_or_derived.to_numpy()
    pairs = tree.query_ball_point(p, r=2*np.sin(np.deg2rad(20/60)/2))
    x["nearby_HI_20arcmin_dv500"] = [
        sum(j != i and np.isfinite(v[i]) and np.isfinite(v[j]) and
            abs(v[i]-v[j]) <= 500 for j in ids) for i,ids in enumerate(pairs)]
    return x

def run():
    OUT.mkdir(parents=True, exist_ok=True)
    receipts = {"start_utc":utc(), "config_sha256":sha((BASE/"freeze.json").read_bytes()),
                "footprints":[], "catalogues":[], "errors":[]}
    mocs=[]
    for i,f in enumerate(CFG["footprint"]):
        raw, rec=fetch_bytes("footprint_"+f["id"], url=f["url"])
        if raw is not None:
            try:
                moc=MOC.from_fits(str(CACHE/("footprint_"+f["id"]+".fits")))
                mocs.append(moc)
                rec["sky_fraction"] = float(moc.sky_fraction)
                rec["max_order"] = int(moc.max_order)
            except Exception as exc:
                rec["status"]="MOC_PARSE_ERROR";rec["error"]=str(exc)[:350]
        receipts["footprints"].append(rec)
    if len(mocs)!=len(CFG["footprint"]):
        receipts["errors"].append("Incomplete official Euclid DR1 footprint; cannot certify overlap")
    allcat=[]
    for cat in CFG["catalogs"]:
        df,rec=input_catalog(cat)
        if df is not None:
            try:
                normed,met=normalize_catalog(df,cat)
                rec.update(met)
                allcat.append(normed)
            except Exception as exc:
                rec["status"]="NORMALIZATION_ERROR";rec["error"]=str(exc)[:600]
        receipts["catalogues"].append(rec)
    del df
    if allcat and len(mocs)==len(CFG["footprint"]):
        # A successful but empty public catalogue response is not a catalogue
        # completion; WALLABY PDR2 has published thousands of detections.
        empty = [c["id"] for c in CFG["catalogs"] if not any(
            len(x)>0 and x.survey.iloc[0]==c["id"] for x in allcat)]
        if empty: receipts["errors"].append("Missing or zero-row surveys: "+",".join(empty))

        full = pd.concat(allcat, ignore_index=True)
        full=coverage(full,mocs)
        match=env_counts(full)
        match=match.sort_values(["survey","id"],kind="stable")
        match.to_csv(OUT/"footprint_hits.csv", index=False, float_format="%.7g")
        # No catalog match is accepted as a stellar mass nondetection.
        match["galaxy_like_screen"] = (
            match.w50_kms.ge(30) & match.vopt_kms_or_derived.ge(600) &
            (match.logMHI.ge(7) | match.logMHI.isna()))
        priority=match[match.galaxy_like_screen].copy()
        priority=priority.sort_values(["nearby_HI_20arcmin_dv500","logMHI"],
                                      ascending=[True,False], na_position="last").head(CFG["followup_cap"])
        priority.to_csv(OUT/"priority_for_followup.csv",index=False,float_format="%.7g")
        counts=[]
        for c in CFG["catalogs"]:
            a=full[full.survey==c["id"]]
            b=match[match.survey==c["id"]]
            counts.append({"survey":c["id"],"input_rows":len(a),
                "valid_positions":int((np.isfinite(a.ra_deg)&np.isfinite(a.dec_deg)).sum()),
                "euclid_center_hits":int((b.euclid_coverage_code==2).sum()),
                "euclid_beam_edge_samples":int((b.euclid_coverage_code==1).sum()),
                "galaxy_like_covered":int((b.galaxy_like_screen).sum()) if "galaxy_like_screen" in b else
                    int((b.w50_kms.ge(30)&b.vopt_kms_or_derived.ge(600)&(b.logMHI.ge(7)|b.logMHI.isna())).sum()),
                "catalog_loaded":True})
        summary={"start_utc":receipts["start_utc"],"complete_all_catalogues":len(allcat)==len(CFG["catalogs"]) and not receipts["errors"],
                 "complete_footprint":True,"counts":counts, "input_total":len(full),
                 "in_footprint_center_or_sampled_edge":len(match),
                 "shortlist_count":len(priority),"stellar_mass_upper_limits_count":0,
                 "confirmed_starless_HI_galaxies":0,
                 "discovery_validation":"NOT_PERFORMED; no physical classification without optical/NIR pixels, deep limits and source vetting",
                 "coverage_warning":"Euclid tile MOC center/16 beam-edge probes, not VIS/NISP pixel coverage or complete beam-polygon intersection"}
    else:
        summary={"complete_all_catalogues":False,"complete_footprint":len(mocs)==len(CFG["footprint"]),
                 "input_total":sum(len(x) for x in allcat),"counts":[],
                 "confirmed_starless_HI_galaxies":0,"stellar_mass_upper_limits_count":0,
                 "discovery_validation":"BLOCKED: required inputs or Euclid MOCs unavailable"}
    receipts["finish_utc"]=utc()
    (OUT/"receipts.json").write_text(json.dumps(receipts,indent=2,allow_nan=False)+"\n")
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(json.dumps(summary,indent=2),flush=True)
    if not summary.get("complete_all_catalogues") or not summary.get("complete_footprint"):
        sys.exit(2)

if __name__=="__main__":
    run()
