#!/usr/bin/env python3
"""Independently recompute survey aggregates from committed bounded batch receipts."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from collections import Counter

import scan

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
BATCHES=OUT/"batches"

def js(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n")

def finite(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None

def snr(r, band):
    a=finite(r.get("flux_"+band+"_3fwhm_aper"))
    b=finite(r.get("fluxerr_"+band+"_3fwhm_aper"))
    return a/b if a is not None and b and b>0 else None

def candidate_rank(r):
    y=finite(r.get("Y_minus_J_AB"))
    jh=finite(r.get("J_minus_H_AB"))
    j=snr(r, "j")
    if y is None:
        y=finite(r.get("Y_minus_J_2sig_limit"))
    if y is None or jh is None or j is None:
        return -1e6
    return round(min(j,25)/5 + min(max(y,0),2) -
                 2*abs(jh-.12) + (2 if r.get("morphology_diagnostic") else -4), 5)

def optical_veto(r):
    # Conservative: a ground-based detection, not proof of low redshift
    # without checking source association and band coverage.
    details={}
    for b in scan.EXT_BANDS:
        q=snr(r,b)
        details[b]=round(q,3) if q is not None else None
    return any(v is not None and v>3 for v in details.values()), details

def main():
    items=scan.catalogue_index()
    batches=math.ceil(len(items)/scan.BATCH_SIZE)
    indices=[]
    for n in range(batches):
        manifest=BATCHES/f"batch_{n:03d}_manifest.json"
        if manifest.exists():
            indices.append(n)
        else:
            break
    # Fail closed on a hole in completed batch prefix.
    if set(indices)!=set(n for n in range(batches)
                         if (BATCHES/f"batch_{n:03d}_manifest.json").exists()):
        raise RuntimeError("Missing intermediate batch")
    total=Counter()
    control_matches=[]
    allrows=[]
    seen_partitions=set()
    for n in indices:
        m=json.loads((BATCHES/f"batch_{n:03d}_manifest.json").read_text())
        p=BATCHES/f"batch_{n:03d}_candidates.jsonl.gz"
        if scan.sha(p.read_bytes())!=m["candidates_sha256"]:
            raise RuntimeError(f"Corrupt candidate receipt batch {n}")
        if m["index_sha256"]!=scan.FROZEN_INDEX_SHA256:
            raise RuntimeError("Mixed publisher index versions")
        batch_count=0
        with gzip.open(p, "rt", encoding="utf-8") as f:
            for line in f:
                row=json.loads(line)
                allrows.append(row)
                batch_count+=1
        if batch_count != m["candidate_count"]:
            raise RuntimeError(f"Candidate count differs from frozen output {n}")
        for part in m["partitions"]:
            key=(part["Norder"],part["Npix"])
            if key in seen_partitions:
                raise RuntimeError("Partition scanned twice: "+str(key))
            seen_partitions.add(key)
            total["rows"]+=part["rows"]
            total["broad"]+=part["broad"]
            total["strict"]+=part["strict"]
            total["pathological_blue"]+=part["blue_JH_below_minus_0p5"]
            control_matches+=part["controls"]
    if total["strict"]!=len(allrows):
        raise RuntimeError("Strict row count mismatch")
    all_ids=[str(r["object_id"]) for r in allrows]
    if len(set(all_ids))!=len(all_ids):
        raise RuntimeError("Duplicate MER object_id in strict selection")
    by_control={}
    for hit in control_matches:
        k=hit["name"]
        if k not in by_control or hit["sep_arcsec"]<by_control[k]["sep_arcsec"]:
            by_control[k]=hit
    controls=[]
    for c in scan.load_controls():
        x=by_control.get(c["name"])
        controls.append({
            "name":c["name"],"class":c["class"],"z":c["z"],
            "nearest":x, "catalog_matched_within_3arcsec":x is not None,
            "strict_recovered":bool(x and x["strict_pass"]),
            "coverage_if_unmatched":"NOT_ESTABLISHED_BY_HEALPIX_PARTITION"})
    matched_q=[c for c in controls if c["class"]=="quasar" and c["nearest"]]
    recovered_q=[c for c in matched_q if c["strict_recovered"]]
    qby_id={str(c["nearest"]["object_id"]):c["name"] for c in matched_q}
    ranked=[]
    n_morph=n_colour=n_optical=0
    for r in allrows:
        known=qby_id.get(str(r["object_id"]))
        jh=finite(r.get("J_minus_H_AB"))
        clean_colour=jh is not None and -.5<=jh<=.45
        clean_morph=bool(r.get("morphology_diagnostic"))
        opt, optical=optical_veto(r)
        n_colour+=int(clean_colour)
        n_morph+=int(clean_morph)
        n_optical+=int(opt)
        if not known and clean_colour and clean_morph and not opt:
            ranked.append({
                "id":str(r["object_id"]), "ra":r["ra"],"dec":r["dec"],
                "Y_minus_J_AB":r["Y_minus_J_AB"],
                "Y_minus_J_2sig_limit":r["Y_minus_J_2sig_limit"],
                "J_minus_H_AB":jh,"J_SNR":round(snr(r,"j"),3),
                "H_SNR":round(snr(r,"h"),3),
                "VIS_SNR":round(snr(r,"vis"),3) if snr(r,"vis") is not None else None,
                "mumax_minus_mag":r.get("mumax_minus_mag"),
                "ellipticity":r.get("ellipticity"),"rank_score":candidate_rank(r),
                "optical_snr":optical,"Norder":r["Norder"],"Npix":r["Npix"],
                "p_quasar_vs_dwarf":None,
                "probability_status":"AWAITING_EMPIRICAL_LABELED_CONTROLS"})
    ranked.sort(key=lambda r:(-r["rank_score"],r["id"]))
    complete=len(indices)==batches
    if complete and (total["rows"]!=scan.PUBLISHED_ROWS or
                     len(seen_partitions)!=len(items)):
        raise RuntimeError("Publisher row/partition coverage validation failed")
    summary={
        "status":"SCAN_COMPLETE" if complete else "PARTIAL_SCAN_NO_SURVEY_CLAIM",
        "index_sha256":scan.FROZEN_INDEX_SHA256,
        "published_wide_rows":scan.PUBLISHED_ROWS,
        "number_of_partitions_in_index":len(items),
        "scanned_partitions":len(seen_partitions),
        "scanned_catalog_rows":total["rows"],"strict_count_in_scanned_partitions":len(allrows),
        "broad_count_in_scanned_partitions":total["broad"],
        "pathologically_blue_JH_below_minus_0p5":total["pathological_blue"],
        "strict_with_inclusive_0p5_lower_JH":n_colour,
        "strict_meeting_nir_mumax_ellipticity_diagnostic":n_morph,
        "strict_with_ground_optical_gt_3sigma":n_optical,
        "ranked_morphology_colour_optical_survivors":len(ranked),
        "controls_yang2026_count":sum(c["class"]=="quasar" for c in controls),
        "controls_yang2026_catalog_matched":len(matched_q),
        "controls_yang2026_strict_recovered":len(recovered_q),
        "recovery_fraction_among_catalog_matched_controls":(
            len(recovered_q)/len(matched_q) if matched_q else None),
        "full_footprint_control_denominator_verified":False,
        "all_published_z_gt_6p5_controls_exhaustively_compiled":False,
        "survey_image_footprint_measured":False,
        "nominal_footprint_deg2":1900,
        "expected_z_gte_7_to_J23_benchmark_approx":19,
        "expected_yang2026_scaled_to_1900deg2_approx":7.6,
        "qlf_expectation_is_completeness_corrected":False,
        "discovery_level":False,
        "limitations":[
            "HATS catalogue objects are not a sky coverage or injection-recovery map",
            "VIS-based POINT_LIKE_PROB is not valid for VIS dropouts",
            "NISP shape criterion transferred from an ultracool-dwarf study is a diagnostic",
            "Strict selection is not a calibrated quasar probability",
            "Yang's 31 discovery controls are not an exhaustive global high-z quasar census",
            "An undetected control may be outside released imaging or missing from the MER catalogue",
            "External spectral and image artefact validation is separate and incomplete"]}
    js(OUT/"summary.json",summary)
    js(OUT/"controls.json",controls)
    with gzip.open(OUT/"ranked_survivors.jsonl.gz", "wt", encoding="utf-8", compresslevel=6) as handle:
        for item in ranked:
            handle.write(json.dumps(item, sort_keys=True, allow_nan=False)+"\n")
    js(OUT/"ranked_candidates.json",ranked[:1000])
    js(OUT/"top_100_unvetted.json",ranked[:100])
    print("AGGREGATE",json.dumps({k:summary[k] for k in
         ("status","scanned_partitions","scanned_catalog_rows",
          "strict_count_in_scanned_partitions",
          "ranked_morphology_colour_optical_survivors",
          "controls_yang2026_catalog_matched","controls_yang2026_strict_recovered")}),
          flush=True)

if __name__=="__main__":
    main()
