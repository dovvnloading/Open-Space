#!/usr/bin/env python3
"""Independent, network-free replay of Euclid HATS receipts and strict-selected rows.

This deliberately does not import scan.py, aggregate.py, or vet.py. It verifies
recorded data rather than copying their implementation decisions. A passed
receipt replay is NOT quasar identification or a footprint completeness test.
"""
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
PART = OUT / "batches"
INDEX_SHA = "f7940bcaed19b6f67295c5592cc81783d798c0f858062d7293e5e3dfbfa57c2a"
EXT = ("g_ext_decam", "r_ext_decam", "i_ext_decam", "z_ext_decam",
       "i_ext_panstarrs", "z_ext_panstarrs")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def val(r, key):
    try:
        x = float(r.get(key))
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def ratio(r, band):
    f = val(r, "flux_" + band + "_3fwhm_aper")
    e = val(r, "fluxerr_" + band + "_3fwhm_aper")
    return f/e if f is not None and e is not None and e > 0 else None


def independent_strict(r):
    v, y, j, h = [val(r, "flux_" + band + "_3fwhm_aper")
                   for band in ("vis", "y", "j", "h")]
    ve, ye, je, he = [val(r, "fluxerr_" + band + "_3fwhm_aper")
                       for band in ("vis", "y", "j", "h")]
    ell = val(r, "ellipticity")
    spur = val(r, "spurious_prob")
    if (None in (v,y,j,h,ve,je,he,ell) or (y<=0 and ye is None) or ve<=0 or je<=0 or he<=0 or
        j<=.65 or h<=0 or not (y<.85*j and j>.65*h) or
        val(r, "vis_det") != 0 or
        (spur is not None and spur>=.7) or
        j/je<6 or h/he<4 or abs(v/ve)>=2 or
        val(r, "flag_j") != 0 or val(r, "flag_h") != 0 or ell>=.3):
        return False
    yj = -2.5*math.log10(y/j) if y>0 else (
        -2.5*math.log10(max(2*ye,y)/j) if max(2*ye,y)>0 else math.inf)
    jh = -2.5*math.log10(j/h)
    return yj>=.4 and jh<=.45


def shape(r):
    m = val(r, "mumax_minus_mag")
    ell = val(r, "ellipticity")
    return m is not None and ell is not None and -3.25<=m<=-2.65 and ell<.2


def optical(r):
    return any((s is not None and s>3) for s in (ratio(r,b) for b in EXT))


def main():
    index_bytes = (OUT/"partition_info.csv").read_bytes()
    assert digest(index_bytes) == INDEX_SHA, "publisher index mismatch"
    index = sorted(csv.DictReader(index_bytes.decode("utf-8-sig").splitlines()),
                   key=lambda x: (int(x["Norder"]),int(x["Npix"])))
    assert len(index)==768 and sum(int(r["num_rows"]) for r in index)==888152165
    summary=json.loads((OUT/"summary.json").read_text())
    assert summary["status"]=="SCAN_COMPLETE", "not full survey"
    assert summary["scanned_partitions"]==768

    stats=Counter()
    ids=set()
    eligible=set()
    control_hits=[]
    scan_source_hash=digest((ROOT/"scan.py").read_bytes())
    for n in range(32):
        path=PART/("batch_%03d_manifest.json" % n)
        manifest=json.loads(path.read_text())
        cells=index[n*24:(n+1)*24]
        assert manifest["batch"]==n and manifest["status"]=="COMPLETE"
        assert manifest["index_sha256"]==INDEX_SHA
        assert manifest["pipeline_sha256"]==scan_source_hash
        assert len(cells)==len(manifest["partitions"])
        compressed=(PART/("batch_%03d_candidates.jsonl.gz" % n)).read_bytes()
        assert digest(compressed)==manifest["candidates_sha256"], "batch hash mismatch"
        for entry,pub in zip(manifest["partitions"],cells):
            assert (entry["Norder"],entry["Npix"]) == (
                int(pub["Norder"]),int(pub["Npix"])), "missing/reordered HATS cell"
            assert entry["rows"]==int(pub["num_rows"]), "publisher row mismatch"
            assert entry["strict"]>=0 and entry["broad"]>=entry["strict"]
            stats["rows"]+=entry["rows"]
            stats["broad"]+=entry["broad"]
            stats["strict_receipts"]+=entry["strict"]
            stats["blue_receipts"]+=entry["blue_JH_below_minus_0p5"]
            control_hits.extend(entry["controls"])
        batch_records=0
        for line in gzip.decompress(compressed).splitlines():
            r=json.loads(line)
            key=str(r["object_id"])
            assert key not in ids, "duplicate strict MER object_id"
            assert independent_strict(r), "a stored strict row fails independent selection"
            assert bool(r["morphology_diagnostic"])==shape(r), "shape flag inconsistent"
            ids.add(key)
            batch_records+=1
            stats["strict_replayed"]+=1
            jh=val(r,"J_minus_H_AB")
            assert jh is not None
            if jh<-.5:
                stats["blue_replayed"]+=1
            if -.5<=jh<=.45:
                stats["inclusive_jh"]+=1
                if shape(r) and not optical(r):
                    eligible.add(key)
            if shape(r):
                stats["shape_rows"]+=1
            if optical(r):
                stats["optical_rows"]+=1
        assert batch_records==manifest["candidate_count"]
    assert stats["rows"]==888152165
    assert stats["strict_replayed"]==stats["strict_receipts"]
    assert stats["blue_replayed"]==stats["blue_receipts"]
    for field,value in (
        ("scanned_catalog_rows",stats["rows"]),
        ("broad_count_in_scanned_partitions",stats["broad"]),
        ("strict_count_in_scanned_partitions",stats["strict_replayed"]),
        ("pathologically_blue_JH_below_minus_0p5",stats["blue_replayed"]),
        ("strict_with_inclusive_0p5_lower_JH",stats["inclusive_jh"]),
        ("strict_meeting_nir_mumax_ellipticity_diagnostic",stats["shape_rows"]),
        ("strict_with_ground_optical_gt_3sigma",stats["optical_rows"])):
        assert summary[field]==value, "summary mismatch: "+field

    nearest={}
    for hit in control_hits:
        name=hit["name"]
        if name not in nearest or hit["sep_arcsec"]<nearest[name]["sep_arcsec"]:
            nearest[name]=hit
    controls=json.loads((OUT/"controls.json").read_text())
    qmatched=[]; qmorph=[]; qstrict=[]; qstrictmorph=[]
    qknown_ids=set()
    for c in controls:
        hit=nearest.get(c["name"])
        assert c["nearest"]==hit
        assert c["catalog_matched_within_3arcsec"]==(hit is not None)
        assert c["strict_recovered"]==bool(hit and hit["strict_pass"])
        if c["class"]=="quasar" and hit:
            qmatched.append(c)
            qknown_ids.add(str(hit["object_id"]))
            m=hit.get("mumax_minus_mag")
            ell=hit.get("ellipticity")
            morph=(m is not None and ell is not None and -3.25<=m<=-2.65 and ell<.2)
            if morph: qmorph.append(c)
            if hit["strict_pass"]:
                qstrict.append(c)
                if morph: qstrictmorph.append(c)
    assert len(controls)==139, "frozen controls changed"
    assert len(qmatched)==summary["controls_combined_catalog_matched"]
    assert len(qstrict)==summary["controls_combined_strict_recovered"]
    ranked=[]
    with gzip.open(OUT/"ranked_survivors.jsonl.gz","rt") as f:
        for line in f: ranked.append(json.loads(line))
    ranked_ids=[str(r["id"]) for r in ranked]
    assert len(ranked_ids)==len(set(ranked_ids))
    assert set(ranked_ids)==eligible-qknown_ids, "ranked candidate cohort differs from strict rows"
    assert len(ranked)==summary["ranked_morphology_colour_optical_survivors"]
    assert ranked==sorted(ranked,key=lambda r:(-r["rank_score"],r["id"]))

    external_file=OUT/"external_receipts.json"
    external=json.loads(external_file.read_text()) if external_file.exists() else []
    model_file=OUT/"probability_model.json"
    model=json.loads(model_file.read_text()) if model_file.exists() else {}
    report={
      "status":"INDEPENDENT_RECEIPT_REPLAY_PASS_NO_DISCOVERY_CLAIM",
      "index_sha256":INDEX_SHA,
      "full_index_partitions_verified":len(index),
      "hats_healpix_cell_union_deg2_not_imaging_area":round(sum(
          41252.96124941927/(12*(4**int(item["Norder"]))) for item in index),3),
      "nominal_DR1_wide_release_area_deg2_about":1900,
      "effective_clean_imaging_area_measured":False,
      "publisher_rows_verified":stats["rows"],
      "strict_rows_independently_checked":stats["strict_replayed"],
      "pathologically_blue_JH_below_minus_0p5":stats["blue_replayed"],
      "ranked_diagnostic_survivors_verified":len(ranked),
      "frozen_quasar_controls":sum(c["class"]=="quasar" for c in controls),
      "catalog_matched_quasar_controls":len(qmatched),
      "strict_recovered_quasar_controls":len(qstrict),
      "strict_recovery_fraction_among_catalog_matched":len(qstrict)/len(qmatched) if qmatched else None,
      "catalog_matched_z_gt_7_controls":sum(c["z"] is not None and c["z"]>7 for c in qmatched),
      "strict_recovered_z_gt_7_controls":sum(c["z"] is not None and c["z"]>7 for c in qstrict),
      "matched_controls_in_borrowed_NISP_shape_band":len(qmorph),
      "strict_recovered_controls_in_borrowed_NISP_shape_band":len(qstrictmorph),
      "external_candidate_receipt_count":len(external),
      "colour_probability_model_status":model.get("model_status","UNAVAILABLE"),
      "catalogue_release_acceptable_if_all_qualifiers_retained":True,
      "discovery_ready":False,
      "referee_limitations":[
        "NISP borrowed dwarf-shape diagnostic excludes nearly all matched known quasars: not a validated quasar-selection completeness criterion",
        "Strict cut itself has poor and nonrepresentative known-quasar recovery",
        "HATS cells and catalogue rows are not a measured usable image area or injection-recovery map",
        "Frozen known-object controls are dated and incomplete for November 2026",
        "Two-class colour posteriors omit galaxies, blending and image artefacts; cannot describe field purity",
        "unWISE/CatWISE, Gaia and literature service errors are not nondetections",
        "No top-object VIS/NISP detector-pixel artefact clearance or new-quasar spectrum"
      ]
    }
    (OUT/"referee.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("INDEPENDENT_REFEREE",json.dumps({
        "pass":True,"partitions":768,"strict":stats["strict_replayed"],
        "controls":len(qstrict),"of":len(qmatched),"morph_controls":len(qmorph),
        "external_receipts":len(external)}),flush=True)


if __name__=="__main__":
    main()
