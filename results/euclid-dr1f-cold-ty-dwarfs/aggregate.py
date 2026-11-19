#!/usr/bin/env python3
"""Independent per-batch accounting, dated control recovery and finite shortlist."""
import gzip
import hashlib
import heapq
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
BATCH=OUT/"batches"
EXPECTED_ROWS=888152165
EXPECTED_PARTITIONS=768
INDEX_SHA="f7940bcaed19b6f67295c5592cc81783d798c0f858062d7293e5e3dfbfa57c2a"
OPTICAL=("g_ext_decam","r_ext_decam","i_ext_decam","z_ext_decam",
         "i_ext_panstarrs","z_ext_panstarrs")

def load(path):
    return json.loads(path.read_text())

def write(name,data):
    (OUT/name).write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+"\n")

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def good(value):
    try:
        value=float(value)
        return value if math.isfinite(value) else None
    except (ValueError,TypeError):return None

def snr(row,band):
    f=good(row.get("flux_"+band+"_3fwhm_aper"))
    e=good(row.get("fluxerr_"+band+"_3fwhm_aper"))
    return f/e if f is not None and e and e>0 else None

def opt(row):
    values={k:snr(row,k) for k in OPTICAL}
    return [k for k,v in values.items() if v is not None and v>3]

def xyz(a,d):
    a,d=np.deg2rad(a),np.deg2rad(d)
    return np.array([np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)]).T

def angular(a,d,b,e):
    return math.hypot(((a-b+180)%360-180)*math.cos(math.radians((d+e)/2)),d-e)*3600

def literature_index():
    with gzip.open(OUT/"known_population.json.gz","rt") as f: literature=json.load(f)
    goodrows=[r for r in literature if good(r.get("ra")) is not None and good(r.get("dec")) is not None]
    tree=cKDTree(xyz([r["ra"] for r in goodrows],[r["dec"] for r in goodrows]))
    return goodrows,tree

def known_matches(r,refs,tree):
    a,d=r["ra"],r["dec"]
    idx=tree.query_ball_point(xyz(a,d),2*math.sin(math.radians(16/3600)/2))
    exact=[];uncertain=[]
    for i in idx:
        item=refs[i]
        sep=angular(a,d,item["ra"],item["dec"])
        is_epoch_robust=bool(item.get("motion_propagated")) or (
            item["source"]=="Kirkpatrick2024_20pc" and item.get("reference_epoch") is not None
            and abs(item["reference_epoch"]-2024.25)<3)
        rec={"catalog":item["source"],"name":item["name"],
             "spectral_type":item["spectral_type"],"sep_arcsec":round(sep,3),
             "epoch_robust":is_epoch_robust}
        if sep<=3.0 and is_epoch_robust:exact.append(rec)
        elif sep<=16:uncertain.append(rec)
    exact.sort(key=lambda x:x["sep_arcsec"])
    uncertain.sort(key=lambda x:x["sep_arcsec"])
    return exact[:5],uncertain[:5]

def score(r):
    y=good(r.get("Y_minus_J_for_selection"))
    h=good(r.get("J_minus_H_for_selection"))
    j=snr(r,"j")
    if y is None or h is None or j is None:return -1000
    # Keep extreme outliers visible but never *promote* arbitrarily blue
    # aperture artefacts as colder; no calibrated posterior.
    ypart=max(0,min(1.6,y-.35))
    hpart=max(0,min(1.2,.15-h))
    extreme=max(0,y-2.2)*2 + max(0,-1.4-h)*2
    psf=good(r.get("mumax_minus_mag"))
    shape=0.2 if psf is not None and -3.25<=psf<=-2.65 else 0
    return round(ypart+1.25*hpart+min(j,30)*.032+shape-extreme,5)

def main():
    indices=[]
    for n in range(32):
        if (BATCH/f"batch_{n:03d}_manifest.json").exists():indices.append(n)
        else:break
    if set(indices)!={int(f.name[6:9]) for f in BATCH.glob("batch_*_manifest.json")}:
        raise RuntimeError("Hole in ordered checkpoints")
    controls=load(OUT/"control_catalog.json")
    control_hits={}
    refs,tree=literature_index()
    totals=Counter()
    partitions=set()
    selected_ids=set()
    heap=[]
    diagnostics=Counter()
    for n in indices:
        m=load(BATCH/f"batch_{n:03d}_manifest.json")
        source=BATCH/f"batch_{n:03d}_candidates.jsonl.gz"
        if m["index_sha256"]!=INDEX_SHA or m["status"]!="COMPLETE":
            raise ValueError("Wrong index or incomplete batch")
        if sha(source)!=m["candidates_sha256"]:
            raise ValueError("Corrupt candidate packet "+str(n))
        if not m.get("selection_sha256"):
            raise ValueError("Missing independent selection fingerprint")
        total_for_batch=0
        with gzip.open(source,"rt",encoding="utf8") as inp:
            for line in inp:
                r=json.loads(line)
                total_for_batch+=1
                ident=str(r["object_id"])
                if ident in selected_ids:
                    raise RuntimeError("Duplicate MER object ID "+ident)
                selected_ids.add(ident)
                h=good(r.get("J_minus_H_for_selection"))
                y=good(r.get("Y_minus_J_for_selection"))
                diagnostics["H_is_limit"]+=bool(r.get("J_minus_H_is_upper_limit"))
                diagnostics["Y_is_limit"]+=bool(r.get("Y_minus_J_is_lower_limit"))
                diagnostics["very_blue_JH_below_minus_1p4"]+=bool(h is not None and h< -1.4)
                veto=opt(r)
                diagnostics["optical_gt_3sigma"]+=bool(veto)
                if y is None or h is None or y<.5 or h>-.15 or veto:continue
                diagnostics["cold_colour_and_optical_survivors"]+=1
                sc=score(r)
                # 250 highest-quality *screening leads*, not 250 discoveries.
                key=(sc,ident)
                new=(sc,ident,{k:r.get(k) for k in (
                    "object_id","ra","dec","Y_minus_J_AB","J_minus_H_AB",
                    "Y_minus_J_for_selection","J_minus_H_for_selection",
                    "Y_minus_J_is_lower_limit","J_minus_H_is_upper_limit",
                    "flux_j_3fwhm_aper","fluxerr_j_3fwhm_aper",
                    "flux_h_3fwhm_aper","fluxerr_h_3fwhm_aper",
                    "flux_y_3fwhm_aper","fluxerr_y_3fwhm_aper",
                    "flux_vis_3fwhm_aper","fluxerr_vis_3fwhm_aper",
                    "mumax_minus_mag","ellipticity","spurious_prob","Norder","Npix")})
                if len(heap)<250:heapq.heappush(heap,new)
                elif key>(heap[0][0],heap[0][1]):heapq.heapreplace(heap,new)
        if total_for_batch!=m["candidate_count"]:
            raise ValueError("Selected-row count mismatch batch "+str(n))
        totals["selected"]+=total_for_batch
        for p in m["partitions"]:
            key=(p["Norder"],p["Npix"])
            if key in partitions:raise RuntimeError("Duplicate partition "+str(key))
            partitions.add(key)
            totals["rows"]+=p["rows"]
            totals["broad"]+=p["broad"]
            totals["strict"]+=p["strict"]
            for hit in p["controls"]:
                name=hit["name"]
                if name not in control_hits or hit["sep_arcsec"]<control_hits[name]["sep_arcsec"]:
                    control_hits[name]=hit
    if totals["strict"]!=totals["selected"]:
        raise RuntimeError("Strict-mask/batch-register mismatch")
    control_details=[]
    for c in controls:
        hit=control_hits.get(c["name"])
        control_details.append({
            "name":c["name"],"class":c["class"],"source":c["source"],
            "spectral_type":c.get("spectral_type"),"predicted_ra_deg":c["ra"],
            "predicted_dec_deg":c["dec"],"motion_propagated":c.get("motion_propagated"),
            "catalog_matched_3arcsec":bool(hit),"mask_recovered":bool(hit and hit["strict_pass"]),
            "nearest":hit,"coverage_verified":False})
    counters={}
    for kind in ("cold_dwarf","warm_T_control","quasar"):
        group=[c for c in control_details if c["class"]==kind]
        matched=[x for x in group if x["catalog_matched_3arcsec"]]
        counters[kind]={"frozen_controls":len(group),"matched_catalog_within_3arcsec":len(matched),
            "mask_recovered":sum(x["mask_recovered"] for x in group),
            "conditional_recovery_if_matched":round(sum(x["mask_recovered"] for x in matched)/len(matched),4) if matched else None,
            "unmatched_coverage_not_verified":len(group)-len(matched)}
    ranked=[]
    for sc,ident,r in sorted(heap,key=lambda x:(-x[0],x[1])):
        exact,uncertain=known_matches(r,refs,tree)
        r.update({"id":str(ident),"score_heuristic":sc,"J_signal_to_noise":round(snr(r,"j"),2),
                  "J_AB":round(23.9-2.5*math.log10(r["flux_j_3fwhm_aper"]),3),
                  "literature_near_epoch_matches":exact,
                  "literature_epoch_uncertain_alerts":uncertain,
                  "known_status":"KNOWN_ASSOCIATION" if exact else (
                      "POSSIBLE_PUBLISHED_ASSOCIATION" if uncertain else "NO_MATCH_IN_DATED_THREE_REFERENCE_SETS"),
                  "not_a_confirmed_cold_dwarf":True})
        ranked.append(r)
    complete=(len(indices)==32)
    if complete and (totals["rows"]!=EXPECTED_ROWS or len(partitions)!=EXPECTED_PARTITIONS):
        raise RuntimeError("Publisher denominator or partition coverage differs")
    summary={"status":"SCAN_COMPLETE" if complete else "PARTIAL_NO_FOOTPRINT_WIDE_CLAIM",
        "data_release":"Euclid DR1-Foundation Wide MER 2026-11-12",
        "source_index_sha256":INDEX_SHA,"advertised_partitions":EXPECTED_PARTITIONS,
        "scanned_partitions":len(partitions),"advertised_rows":EXPECTED_ROWS,
        "scanned_rows":totals["rows"],"broad_mask_rows":totals["broad"],
        "colour_mask_rows":totals["strict"],"strict_rows_in_checkpoint":totals["selected"],
        "diagnostics":dict(diagnostics),
        "control_recovery":counters,
        "pre_external_shortlist_length":len(ranked),
        "pre_external_shortlist_without_exact_literature_match":sum(not r["literature_near_epoch_matches"] for r in ranked),
        "not_independently_measured":{"unmasked_area":True,"NISP_PSF_star_galaxy":True,
                                     "Euclid_image_pixel_artifacts":True,
                                     "photometric_injection_recovery":True},
        "discovery_level":False}
    write("summary.json",summary)
    write("controls_recovery.json",control_details)
    write("unvetted_shortlist.json",ranked)
    print("AGGREGATE",json.dumps({k:summary[k] for k in ("status","scanned_partitions",
          "scanned_rows","colour_mask_rows","pre_external_shortlist_length","control_recovery")}),flush=True)

if __name__=="__main__":
    main()
