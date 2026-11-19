#!/usr/bin/env python3
"""Frozen eRASS1 DR1 catalogue screen. Network and computations run ONLY on Actions."""
from __future__ import annotations
import csv
import hashlib
import io
import json
import math
import os
import sys
import tarfile
import tempfile
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
OUT.mkdir(parents=True, exist_ok=True)
MAIN = "https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Main.v1.2.fits.tar.gz"
HARD = "https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Hard.tar.gz"
XMATCH = "https://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"
# Freeze these before seeing catalogue ranks. No threshold tuning to controls.
PROTOCOL = {
    "release": "eROSITA-DE eRASS1 DR1 main v1.2 (2026-01-12); hard official DR1",
    "main_expected_rows": 930203, "hard_expected_rows": 5466,
    "min_det_like": 12.0, "min_ml_counts": 20.0, "max_poserr_arcsec": 10.0,
    "max_extent_arcsec": 0.0, "max_ext_like": 6.0,
    "min_soft_fraction_P1P2_over_P1P2P3": 0.90,
    "min_P1_counts": 5.0, "min_flux_0p2_2p3": 1e-13,
    "exclude_flags": ["FLAG_OPT", "FLAG_SP_BPS", "FLAG_SP_SCL", "FLAG_SP_LGA", "FLAG_SP_GC_CONS"],
    "exclude_strong_hard_association": True,
    "nearest_brighter_neighbor_arcsec": 60.0, "brighter_neighbor_factor": 5.0,
    "crossmatch_radius_arcsec": 8.0, "crossmatch_top_n": 60,
    "raw_screen_output_cap": 20000,
    "optic_to_X_ratio_required": 1e4,
    "ratio_definition": "Maccacaro log10(f_X/f_V) = log10(f_X[erg cm-2 s-1]) + 0.4 V + 5.37; g=V is an UNVALIDATED proxy",
    "eligibility": "strict heuristic; not a completeness-calibrated all-sky selection",
}
SURVEYS = {
    "Gaia_DR3": "I/355/gaiadr3",
    "PanSTARRS_DR1": "II/349/ps1",
    "DES_DR2": "II/371/des_dr2",
    "unWISE": "II/363/unwise",
    "CatWISE2020": "II/365/catwise",
    "2MASS": "II/246/out",
    "XMM_4XMM": "IX/59/xmm4dr9",
    "Chandra_CSC2": "IX/57/csc2",
    "Swift_2SXPS": "IX/58/2sxpscle",
    "ATNF_VizieR_2016": "B/psr/psr",
}
OPTICAL = list(SURVEYS)[:6]
ARCHIVES = ["XMM_4XMM", "Chandra_CSC2", "Swift_2SXPS"]
# Reference positions are for POSITIVE CONTROLS only, not search seeds.
# Literature positions are epoch-dependent for high-proper-motion stars.
CONTROL_SEXAGESIMAL = {
    "RX_J0420.0-5022": ("04 20 01.94", "-50 22 47.8"),
    "RX_J0720.4-3125": ("07 20 24.93", "-31 25 49.78"),
    "RX_J0806.4-4123": ("08 06 23.40", "-41 22 30.9"),
    "RX_J1308.6+2127_RBS1223": ("13 08 48.27", "+21 27 06.8"),
    "RX_J1605.3+3249": ("16 05 18.52", "+32 49 18.0"),
    "RX_J1856.5-3754": ("18 56 35.795", "-37 54 35.54"),
    "RX_J2143.0+0654_RBS1774": ("21 43 03.4", "+06 54 17.5"),
    "eRASSU_J131716.9-402647_published": ("13 17 16.9", "-40 26 47"),
    "eRASSU_J065715.3+260428_published": ("06 57 15.3", "+26 04 28"),
}
RECEIPTS = {"created_utc": datetime.now(timezone.utc).isoformat(), "protocol": PROTOCOL,
            "requests": [], "catalogues": [], "literature": [], "status": "started"}
STEPS = {"errors": [], "completed": [], "pending": [], "counts": {}, "claims": []}
def save(name, value):
    (OUT/name).write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8*1024*1024), b""):
            h.update(block)
    return h.hexdigest()

def download(session, label, url, dest):
    record = {"label": label, "url": url, "utc": datetime.now(timezone.utc).isoformat()}
    for attempt in range(1, 4):
        try:
            with session.get(url, stream=True, timeout=(30, 180)) as r:
                record.update({"http_status": r.status_code, "final_url": r.url,
                               "content_type": r.headers.get("Content-Type"),
                               "last_modified": r.headers.get("Last-Modified"),
                               "etag": r.headers.get("ETag"), "attempt": attempt})
                r.raise_for_status()
                with open(dest, "wb") as f:
                    for chunk in r.iter_content(4*1024*1024):
                        if chunk: f.write(chunk)
            record["bytes"] = dest.stat().st_size
            record["sha256"] = sha256(dest)
            RECEIPTS["requests"].append(record)
            return
        except Exception as e:
            record["error"] = repr(e)
            if attempt < 3: time.sleep(5 * attempt)
    RECEIPTS["requests"].append(record)
    raise RuntimeError(label + " download failed: " + str(record))

def unpack_one(archive, directory, label):
    with tarfile.open(archive, "r:*") as t:
        members = [m for m in t.getmembers() if m.isfile() and m.name.lower().endswith((".fits", ".fit", ".fits.gz"))]
        if len(members) != 1:
            raise RuntimeError(label + ": expected one FITS member; got " + repr([(m.name, m.size) for m in members]))
        m = members[0]
        # Do not trust paths in remote tar files.
        target = directory/(label + (".fits.gz" if m.name.endswith(".gz") else ".fits"))
        with t.extractfile(m) as src, open(target, "wb") as dest:
            while True:
                chunk = src.read(8*1024*1024)
                if not chunk: break
                dest.write(chunk)
        return target, {"member": m.name, "fits_bytes": target.stat().st_size, "fits_sha256": sha256(target)}

def num(data, name, default=0.):
    if name not in data.names:
        raise RuntimeError("Required official main catalogue field absent: " + name)
    x = np.asarray(data[name], dtype=np.float64)
    return np.nan_to_num(x, nan=default, posinf=default, neginf=default)

def vec(coords):
    return np.array(coords.cartesian.xyz.value).T

def write_csv(path, rows, fields):
    with open(OUT/path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows: writer.writerow({k: row.get(k, "") for k in fields})

def references(main_data, radec, selected_indices, coords):
    results = {}
    all_coords = SkyCoord(radec[0]*u.deg, radec[1]*u.deg)
    known = SkyCoord([p[0] for p in CONTROL_SEXAGESIMAL.values()],
                     [p[1] for p in CONTROL_SEXAGESIMAL.values()], unit=(u.hourangle, u.deg), frame="icrs")
    ix, sep, _ = known.match_to_catalog_sky(all_coords)
    selected = set(int(x) for x in selected_indices)
    for (name, _), k, d, point in zip(CONTROL_SEXAGESIMAL.items(), ix, sep, known):
        i = int(k)
        results[name] = {
            "reference_ra": round(point.ra.deg, 6), "reference_dec": round(point.dec.deg, 6),
            "reference_galactic_l": round(float(point.galactic.l.deg), 3),
            "rough_western_hemisphere": bool(point.galactic.l.deg > 180),
            "nearest_IAUNAME": str(main_data["IAUNAME"][i]).strip(),
            "nearest_catalogue_row": i,
            "separation_arcsec": round(float(d.arcsec), 3),
            "detected_within_60arcsec": bool(d.arcsec <= 60),
            "passed_frozen_screen": bool(i in selected),
            "proper_motion_and_epoch_warning": True,
        }
    return results

def xmatch_catalog(session, name, catalogue, sources):
    if not sources: return {"status":"not_run_empty"}, {}
    payload = io.StringIO()
    w = csv.DictWriter(payload, fieldnames=["srcid", "ra", "dec"])
    w.writeheader()
    for row in sources: w.writerow({"srcid": row["uid"], "ra": row["ra"], "dec": row["dec"]})
    info = {"catalogue": catalogue, "service": XMATCH, "source_count": len(sources),
            "max_arcsec": PROTOCOL["crossmatch_radius_arcsec"], "requested_utc": datetime.now(timezone.utc).isoformat()}
    try:
        with session.post(XMATCH, data={
            "request": "xmatch", "distMaxArcsec": str(PROTOCOL["crossmatch_radius_arcsec"]),
            "RESPONSEFORMAT": "csv", "cat2": "vizier:" + catalogue,
            "colRA1": "ra", "colDec1": "dec"
        }, files={"cat1": ("src.csv", payload.getvalue().encode(), "text/csv")}, timeout=(30, 160)) as r:
            info["http_status"] = r.status_code
            info["bytes"] = len(r.content)
            info["response_sha256"] = hashlib.sha256(r.content).hexdigest()
            info["final_url"] = r.url
            r.raise_for_status()
            lines = r.text
        if lines.lstrip().startswith(("ERROR", "<", "{")):
            raise RuntimeError("service did not return CSV: " + lines[:220])
        table = list(csv.DictReader(io.StringIO(lines)))
        match = {}
        for row in table:
            try:
                key = str(int(float(row["srcid"])))
                arc = float(row.get("angDist", "nan"))
            except (KeyError, ValueError, TypeError): continue
            if not math.isfinite(arc): continue
            if key not in match or arc < match[key]["separation_arcsec"]:
                match[key] = {"separation_arcsec": round(arc,4),
                              "catalogue_id": next((str(row[k])[:70] for k in ["Source","_2MASS","objID","Name","IAUName","PS1"] if row.get(k)), ""),
                              "match_row_sha256": hashlib.sha256(json.dumps(row,sort_keys=True).encode()).hexdigest()}
        info.update({"status": "ok", "returned_rows": len(table), "matched_sources": len(match)})
        return info, match
    except Exception as e:
        info.update({"status":"failed", "error": repr(e)})
        return info, {}

def match_published_kurpas():
    # The VizieR paper is independent of the fresh search. A failed query is NOT a novelty pass.
    info = {"source": "Kurpas et al. 2024 A&A 687 A251",
            "url": "https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/687/A251"}
    try:
        from astroquery.vizier import Vizier
        vz = Vizier(columns=["**"], row_limit=-1, timeout=120)
        tabs = vz.get_catalogs("J/A+A/687/A251")
        info["tables"] = {str(k):len(t) for k,t in tabs.items()}
        rows = []
        for t in tabs.values():
            for row in t:
                vals = {k: str(row[k]).strip() for k in t.colnames}
                ra = next((vals[k] for k in ["RAJ2000","_RAJ2000","RA_ICRS","RAdeg","RA"] if k in vals),None)
                dec = next((vals[k] for k in ["DEJ2000","_DEJ2000","DE_ICRS","DEdeg","Dec","DE"] if k in vals),None)
                if ra is None or dec is None: continue
                try:
                    if ":" in ra or " " in ra: p=SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
                    else: p=SkyCoord(float(ra)*u.deg,float(dec)*u.deg)
                    rows.append({"ra": float(p.ra.deg), "dec":float(p.dec.deg), "source_table":str(t.meta.get("name",""))})
                except (ValueError,TypeError): continue
        info["parsed_positions"] = len(rows)
        info["status"] = "ok" if rows else "no_parseable_positions"
        return info, rows
    except Exception as e:
        info.update({"status": "unavailable", "error":repr(e)})
        return info, []

def run():
    save("protocol.json", PROTOCOL)
    with requests.Session() as session, tempfile.TemporaryDirectory(prefix="erass1-dr1-") as tmp:
        session.headers.update({"User-Agent":"Open-Space-XDINS-research/1.0 (reproducible-academic-screen)"})
        tmp = Path(tmp)
        files = {}
        for name, url, expected in [("main", MAIN,930203), ("hard", HARD,5466)]:
            archive = tmp/(name+".tar.gz")
            download(session, "eRASS1_"+name, url, archive)
            path, meta = unpack_one(archive,tmp,name)
            files[name] = path
            RECEIPTS["catalogues"].append({"name":name,"url":url,
                "archive_sha256":sha256(archive), "archive_bytes":archive.stat().st_size, **meta,
                "expected_rows":expected})
        with fits.open(files["main"], memmap=True, lazy_load_hdus=True) as hdul, fits.open(files["hard"], memmap=True, lazy_load_hdus=True) as hard:
            main = hdul[1].data
            hard_data = hard[1].data
            N = len(main)
            RECEIPTS["catalogues"][0].update({"rows":N, "columns":len(main.names), "hdu_name":hdul[1].name})
            RECEIPTS["catalogues"][1].update({"rows":len(hard_data),"columns":len(hard_data.names),"hdu_name":hard[1].name})
            if N != PROTOCOL["main_expected_rows"] or len(hard_data) != PROTOCOL["hard_expected_rows"]:
                raise RuntimeError("Catalogue version/row-count mismatch; do not classify results")
            STEPS["completed"].append("official eRASS1 main+hard files downloaded, hashed, row counts checked")
            needed = ["RA","DEC","ML_CTS_1","ML_CTS_P1","ML_CTS_P2","ML_CTS_P3",
                      "ML_FLUX_1","POS_ERR","DET_LIKE_0","EXT","EXT_LIKE","IAUNAME","UID","UID_Hard"]
            missing = [k for k in needed if k not in main.names]
            if missing: raise RuntimeError("Expected fields missing: "+repr(missing))
            ra, dec = num(main,"RA"),num(main,"DEC")
            fx=num(main,"ML_FLUX_1"); counts=num(main,"ML_CTS_1")
            p1,p2,p3=[np.maximum(num(main,"ML_CTS_"+k),0) for k in ("P1","P2","P3")]
            softness = (p1+p2)/np.maximum(p1+p2+p3,1e-12)
            flags = np.zeros(N,dtype=bool)
            for name in PROTOCOL["exclude_flags"]:
                if name in main.names: flags |= (num(main,name)!=0)
                else: STEPS["pending"].append("missing upstream flag "+name)
            quality = (np.isfinite(ra)&np.isfinite(dec)&(fx>0))
            det = quality & (num(main,"DET_LIKE_0")>=PROTOCOL["min_det_like"]) & (counts>=PROTOCOL["min_ml_counts"]) & (fx>=PROTOCOL["min_flux_0p2_2p3"])
            shape = det & (num(main,"POS_ERR")>0)&(num(main,"POS_ERR")<=PROTOCOL["max_poserr_arcsec"]) & (num(main,"EXT")<=PROTOCOL["max_extent_arcsec"]) & (num(main,"EXT_LIKE")<=PROTOCOL["max_ext_like"])
            soft = shape & (softness>=PROTOCOL["min_soft_fraction_P1P2_over_P1P2P3"]) & (p1>=PROTOCOL["min_P1_counts"])
            not_flag = soft & ~flags
            # Hard catalogue main association uses sign convention; verify raw hard association counts.
            uid_hard = np.asarray(main["UID_Hard"],dtype=np.int64)
            selected = not_flag & (uid_hard<=0)
            stages={"all_main":int(N),"quality_finite":int(quality.sum()),"flux_det_counts":int(det.sum()),
                    "pointlike":int(shape.sum()),"soft":int(soft.sum()),"upstream_flags_clean":int(not_flag.sum()),
                    "without_strong_hard_match":int(selected.sum()),"hard_rows":len(hard_data)}
            STEPS["counts"].update(stages)
            STEPS["completed"].append("frozen main-catalogue 0.2-2.3keV, P1/P2/P3 and hard-association screen")
            selected_idx = np.flatnonzero(selected)
            controls = references(main,(ra,dec),selected_idx,None)
            save("controls.json", controls)
            # Bright neighbour and confusion flag from all main entries, using spherical 3D geometry.
            eligible = np.flatnonzero(selected)
            all_coord = SkyCoord(ra*u.deg, dec*u.deg)
            all_points = vec(all_coord)
            tree = cKDTree(all_points)
            rows=[]
            for i in eligible:
                near = tree.query_ball_point(all_points[i],r=2*np.sin(np.deg2rad(PROTOCOL["nearest_brighter_neighbor_arcsec"]/3600)/2))
                bright = [j for j in near if j!=i and fx[j]>PROTOCOL["brighter_neighbor_factor"]*fx[i]]
                g_req=(4-math.log10(float(fx[i]))-5.37)/.4
                rows.append({
                    "uid":str(int(main["UID"][i])),
                    "iauname":str(main["IAUNAME"][i]).strip(),
                    "ra":round(float(ra[i]),7),"dec":round(float(dec[i]),7),
                    "flux_0p2_2p3":float(fx[i]),
                    "ml_counts":round(float(counts[i]),3),
                    "det_like":round(float(num(main,"DET_LIKE_0")[i]),3),
                    "soft_fraction":round(float(softness[i]),4),
                    "p1_counts":round(float(p1[i]),2),
                    "poserr_arcsec":round(float(num(main,"POS_ERR")[i]),3),
                    "hard_uid":int(uid_hard[i]), "brighter_neighbors_60arcsec":len(bright),
                    "illustrative_g_equals_V_needed_for_ratio_1e4":round(g_req,2),
                    "catalogue_row":int(i)
                })
            # Fixed ranking by high soft counts and soft fraction, not by optical coincidences.
            rows.sort(key=lambda x:(-x["ml_counts"]*x["soft_fraction"],-x["det_like"],x["uid"]))
            write_csv("screened_main.csv", rows[:PROTOCOL["raw_screen_output_cap"]],
                      ["uid","iauname","ra","dec","flux_0p2_2p3","ml_counts","det_like",
                       "soft_fraction","p1_counts","poserr_arcsec","hard_uid",
                       "brighter_neighbors_60arcsec","illustrative_g_equals_V_needed_for_ratio_1e4","catalogue_row"])
            STEPS["counts"]["screened_frozen_before_neighbor"] = len(rows)
            STEPS["counts"]["screened_output_rows"] = min(len(rows),PROTOCOL["raw_screen_output_cap"])
            clean = [r for r in rows if r["brighter_neighbors_60arcsec"]==0]
            follow = clean[:PROTOCOL["crossmatch_top_n"]]
            STEPS["counts"]["neighbor_clean"] = len(clean)
            STEPS["counts"]["followup_top_n"] = len(follow)
            save("candidate_screen.json", {"stages":stages, "counts":STEPS["counts"],
                "ranking": "descending ML_CTS_1 × soft_fraction; deterministic UID tie-break",
                "followup_capped":len(clean)>len(follow),
                "bright_neighbor": "f_X>5x within 60arcsec in main; flag, not verified image morphology",
                "flux_conversion_warning":"catalogue ML_FLUX_1 assumes catalog spectral model, not a fitted NS blackbody",
                "normalization": "Maccacaro V-band proxy gives a REQUIRED g=V magnitude, NOT a measured limit"})
            STEPS["completed"].append("recorded full frozen-screen denominator and rank-capped follow-up")
        # The complete optical match process is applied to a bounded independent ranking.
        by_name = {}
        for name in [*OPTICAL,*ARCHIVES,"ATNF_VizieR_2016"]:
            receipt, matches = xmatch_catalog(session, name, SURVEYS[name], follow)
            RECEIPTS["requests"].append({"label":"CDS_xmatch_"+name,**receipt})
            by_name[name] = matches
            if receipt["status"]=="ok": STEPS["completed"].append(name+" positional catalogue check")
            else: STEPS["pending"].append(name+" query failed or unverified")
        lrec, published = match_published_kurpas()
        RECEIPTS["literature"].append(lrec)
        if published:
            pc = SkyCoord([r["ra"] for r in published]*u.deg,[r["dec"] for r in published]*u.deg)
            fc = SkyCoord([r["ra"] for r in follow]*u.deg,[r["dec"] for r in follow]*u.deg) if follow else None
            if fc is not None:
                idx, d, _=fc.match_to_catalog_sky(pc)
                for item,i,sep in zip(follow,idx,d):
                    item["kurpas_2024_nearest_arcsec"]=round(float(sep.arcsec),2)
                    item["kurpas_2024_match_30arcsec"]=bool(sep.arcsec<30)
            STEPS["completed"].append("Kurpas 2024 positional candidate comparison where parseable")
        else: STEPS["pending"].append("Kurpas 2024 full list of 33 candidates not independently parsed")
        # LS DR10 query is deliberately a separate service: not interchangeable with Pan-STARRS/DES.
        # VizieR lacks a verified full LS DR10 table here; a missing LS result is never a non-detection.
        STEPS["pending"].append("Legacy Surveys DR10 reliable positional catalog+depth check; photo-z and forced flux")
        STEPS["pending"].extend([
            "survey coverage, epoch and per-position limiting magnitudes for all optical/IR non-matches",
            "source-level coronal-star/AGN/radio/IR-color/cluster classification and non-detection significance",
            "XMM/Chandra/Swift native event-level spectral fits and timing tests; catalogue hits alone are insufficient",
            "current ATNF update (VizieR ATNF table is 2016); SIMBAD and complete 2026 literature novelty vetting",
            "eRASS1 scan-by-scan light curves; crossmatch to Boller et al. 2025 variability catalogue",
            "realistic positive-control and injection/recovery completeness calibration",
        ])
        cfields=["uid","iauname","ra","dec","flux_0p2_2p3","ml_counts","det_like","soft_fraction",
                 "brighter_neighbors_60arcsec","illustrative_g_equals_V_needed_for_ratio_1e4",
                 "optical_surveys_successful","optical_positional_matches","unmatched_with_coverage_unverified",
                 "archival_catalogue_matches","known_pulsar_2016_match",
                 "kurpas_2024_match_30arcsec","status"]
        for item in follow:
            uid=item["uid"]
            succeeded=[name for name in OPTICAL if any(x.get("label")=="CDS_xmatch_"+name and x.get("status")=="ok" for x in RECEIPTS["requests"])]
            matches=[name for name in OPTICAL if uid in by_name[name]]
            archives=[name for name in ARCHIVES if uid in by_name[name]]
            item.update({
                "optical_surveys_successful":";".join(succeeded),
                "optical_positional_matches":";".join(matches),
                "unmatched_with_coverage_unverified":";".join(name for name in succeeded if name not in matches),
                "archival_catalogue_matches":";".join(archives),
                "known_pulsar_2016_match":bool(uid in by_name["ATNF_VizieR_2016"]),
                "status":"KNOWN_KURPAS" if item.get("kurpas_2024_match_30arcsec") else (
                        "COUNTERPART_OR_BLENDING_REVIEW" if matches else "UNCLASSIFIED_COVERAGE_PENDING")
            })
        write_csv("candidates.csv", follow,cfields)
        save("matches.json", {"catalogue_queries":{name:next((x for x in RECEIPTS["requests"] if x.get("label")=="CDS_xmatch_"+name),{}) for name in SURVEYS},
                              "nearest_matches":{name:by_name[name] for name in SURVEYS},
                              "scope":"positional 8arcsec, NOT probabilistic IDs or upper-limit proofs"})
        STEPS["counts"].update({"followup_with_any_optical_match":sum(bool(r["optical_positional_matches"]) for r in follow),
                               "followup_without_any_optical_match":sum(not bool(r["optical_positional_matches"]) for r in follow),
                               "qualified_novel_XDINS":0})
        STEPS["claims"]=["No new XDINS is validated. Do not infer a calibrated null from limited follow-up.",
                          "Zero 'qualified_novel_XDINS' means zero objects cleared all required checks, not no astrophysical objects."]
        STEPS["completed"].append("bounded 8arcsec multi-survey positional checks with per-service HTTP/SHA receipts")
    RECEIPTS["status"]="complete_screen_incomplete_astrophysical_validation"

if __name__=="__main__":
    try:
        run()
    except Exception as exc:
        RECEIPTS["status"]="failed_incomplete"
        STEPS["errors"].append({"type":type(exc).__name__, "message":str(exc), "trace":traceback.format_exc()[-3500:]})
    finally:
        RECEIPTS["finished_utc"]=datetime.now(timezone.utc).isoformat()
        save("receipts.json", RECEIPTS)
        save("checks.json", STEPS)
        print(json.dumps({"receipt_status":RECEIPTS["status"],"counts":STEPS["counts"],"errors":[e["message"] for e in STEPS["errors"]]},sort_keys=True))
    if STEPS["errors"]: sys.exit(2)
