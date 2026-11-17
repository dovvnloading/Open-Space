#!/usr/bin/env python3
"""Independent stdlib-only replay of deposited compact results, NOT native FITS/Parquet."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FROZEN_TILES = "bfd3724847646662e651363d4624728b69cb9a50577548320e5d7fddf89494b6"
FROZEN_RECEIPTS = "a987e71517e3513dabfa6b6b3dae46ccaf1b69cbb3ab6356bb0dd98c557c1e7b"


def rows(path):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def blob(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def total(rs, col):
    return sum(int(r[col]) for r in rs)


def check_coverage():
    tiles = rows("euclid/q1_tiles.csv")
    cutouts = rows("euclid/q1_cutouts.csv")
    cover = rows("detector/tables/coverage.csv")
    summary = blob("detector/tables/summary.json")
    for filename, digest in (("euclid/q1_tiles.csv", FROZEN_TILES),
                             ("euclid/q1_cutouts.csv", FROZEN_RECEIPTS)):
        require(hashlib.sha256((ROOT / filename).read_bytes()).hexdigest() == digest,
                "frozen SHA-256: " + filename)
    require(len(tiles) == len({r["tile"] for r in tiles}) == 352, "352 source tiles")
    require(Counter(r["field"] for r in tiles) ==
            {"EDF-S": 148, "EDF-N": 124, "EDF-F": 72, "LDN1641": 8},
            "full tile inventory by field")
    blind = {r["tile"] for r in tiles if r["field"] != "LDN1641"}
    require(len(blind) == 344, "344 extragalactic centres")
    require(len(cutouts) == 700 and
            Counter(r["product"] for r in cutouts) ==
            {"VIS": 344, "rms": 344, "Y": 4, "J": 4, "H": 4},
            "688 VIS/RMS plus 12 NIR receipts")
    index = {(r["tile"], r["product"]): r for r in cutouts}
    require(len(index) == 700 and
            {r["tile"] for r in cutouts} == blind, "unique receipt keys and tiles")
    for r in cutouts:
        h = r["sha256"]
        require(len(h) == 64 and all(c in "0123456789abcdef" for c in h),
                "bad cutout hash")
        require(int(r["bytes"]) > 0 and int(r["nx"]) > 0 and int(r["ny"]) > 0,
                "bad receipt dimensions/size")
        require(r["tile"] in r["url"], "receipt URL/tile mismatch")
    require(len(cover) == 344 and {r["tile"] for r in cover} == blind,
            "coverage tile set")
    for r in cover:
        v, e = index[(r["tile"], "VIS")], index[(r["tile"], "rms")]
        require((r["vis_sha256"], r["rms_sha256"]) ==
                (v["sha256"], e["sha256"]), "coverage hash/receipt mismatch")
        require((r["nx"], r["ny"]) == (v["nx"], v["ny"]) ==
                (e["nx"], e["ny"]), "science/rms/coverage size mismatch")
    require(sum(int(r["valid_pixels"]) > 0 for r in cover) == 272,
            "272 valid stamp centres")
    require(total(cover, "valid_pixels") == 73152512 and
            total(cover, "all_angles_center_pixels") == 24610032,
            "valid and conservative pixel denominators")
    area = sum(int(r["valid_pixels"]) * float(r["square_arcsec_per_pixel"])
               for r in cover) / (3600**2)
    safe = sum(int(r["all_angles_center_pixels"]) *
               float(r["square_arcsec_per_pixel"]) for r in cover) / (3600**2)
    require(abs(area - summary["valid_area_deg2"]) < 1e-12 and
            abs(safe - summary["full_stripe_center_area_deg2"]) < 1e-12,
            "recomputed footprints differ from summary")
    require(summary["searched_tiles"] == summary["inverse_tiles"] == 344 and
            summary["with_valid_pixels"] == 272 and
            summary["min_center_separation_arcsec"] >
            summary["max_stamp_diagonal_arcsec"], "selection/overlap checks")
    require(tuple(total(cover, x) for x in
                  ("positive_raw", "positive_profile", "negative_raw",
                   "negative_profile")) == (8264, 651, 2529, 450),
            "signed morphology counts")
    for name in ("positive_raw", "positive_profile", "negative_raw",
                 "negative_profile"):
        require(total(cover, name) == summary[name], "summary " + name)
    return area, safe


def check_recovery():
    trials = rows("detector/tables/injections.csv")
    aggregate = rows("detector/tables/recovery.csv")
    by_sample = rows("detector/tables/recovery_by_sample.csv")
    require(len(trials) == 840 and len({r["tile"] for r in trials}) == 27,
            "injection tile/sample counts")
    fields = ("preexisting_match", "detected", "new_match", "new_profile_match")
    bins = defaultdict(list)
    samples = defaultdict(list)
    for r in trials:
        k = (r["length_arcsec"], r["fwhm_arcsec"], r["mu_ab_arcsec2"])
        bins[k].append(r)
        samples[(r["selection"], "all", *k)].append(r)
        if r["strip_valid"] == "True":
            samples[(r["selection"], "valid strip", *k)].append(r)
    expected = {(6.0, 23.0): 62, (6.0, 24.0): 54,
                (6.0, 25.0): 40, (6.0, 26.0): 11, (6.0, 27.0): 0,
                (10.0, 23.0): 64, (10.0, 24.0): 58,
                (10.0, 25.0): 47, (10.0, 26.0): 23, (10.0, 27.0): 2}
    require(len(bins) == len(aggregate) == 10 and len(by_sample) == 40,
            "full injection grid and sample cells")
    for table, lookup, subset in ((aggregate, bins, False),
                                  (by_sample, samples, True)):
        for r in table:
            k = ((r["selection"], r["coverage"]) if subset else ()) + (
                r["length_arcsec"], r["fwhm_arcsec"], r["mu_ab_arcsec2"])
            t = lookup[k]
            require(int(r["attempts"]) == len(t), "trial counts " + str(k))
            for f in fields:
                require(int(r[f]) == sum(x[f] == "True" for x in t),
                        "trial output " + f + " " + str(k))
            if not subset:
                require(len(t) == 84 and float(r["fwhm_arcsec"]) == .2,
                        "84 trials with 0.2 arcsec FWHM per bin")
                require(int(r["new_profile_match"]) ==
                        expected[(float(r["length_arcsec"]),
                                  float(r["mu_ab_arcsec2"]))],
                        "README recovery numerator")
    require(total(aggregate, "new_profile_match") == 361, "injection matches")
    return len(trials)


def check_catalog_and_reviews():
    short = rows("detector/tables/profile_shortlist.csv")
    review = rows("detector/review_rows.csv")
    require(len(short) == 651 and len(review) == 68, "shortlist/review lengths")
    ids = {f'{r["tile"]}_positive_{r["rank"]}' for r in short}
    require(len({r["id"] for r in review}) == 68 and
            {r["id"] for r in review} <= ids, "review IDs trace to shortlist")
    classes = Counter(r["visual_class"] for r in review)
    require(classes["weak or confused"] == 32 and
            classes["stellar diffraction spike"] == 18, "visual classes")
    unresolved = set()
    for line in (ROOT / "detector/unresolved_followup.md").read_text().splitlines():
        if line.startswith("| ") and "positive_" in line:
            item = line.split("|")[1].strip().strip(chr(96))
            require(item not in unresolved, "duplicated unresolved follow-up")
            unresolved.add(item)
    require(unresolved == {r["id"] for r in review
                           if r["visual_class"] == "weak or confused"},
            "unresolved follow-up inventory")
    partitions = rows("catalog/hats_partitions.csv")
    summary = blob("catalog/hats_summary.json")
    selected = rows("catalog/hats_selected.csv")
    gaia = rows("catalog/gaia_alignment.csv")
    require(len(partitions) == summary["expected_partitions"] ==
            summary["partitions"] == 114, "HATS parent partition inventory")
    require(tuple(total(partitions, f) for f in ("n", "n_vis", "n_selected")) ==
            (29953430, 22349951, 5134), "full catalogue denominators")
    require(tuple(summary[f] for f in
                  ("n", "n_vis", "n_selected", "unique_selected_ids")) ==
            (29953430, 22349951, 5134, 5134), "catalogue summary")
    require(len(selected) == len({r["object_id"] for r in selected}) == 5134,
            "selected row count/uniqueness")
    require(all(int(r["mer_vis_det"]) == 1 and
                .95 < float(r["mer_ellipticity"]) < 1 and
                float(r["mer_semimajor_axis"]) >= 7.5 for r in selected),
            "frozen catalogue shape cuts")
    require(len(gaia) == 5134 and {r["object_id"] for r in gaia} ==
            {r["object_id"] for r in selected}, "Gaia join keys")
    require(Counter(r["field"] for r in gaia) ==
            {"EDF-N": 2482, "EDF-S": 1903, "EDF-F": 735, "LDN1641": 14},
            "selected rows per field")
    deep = [r for r in gaia if r["field"] != "LDN1641"]
    require(sum(r["aligned"] == "True" for r in deep) == 4490 and
            sum(r["rotated"] == "True" for r in deep) == 92,
            "deep-field Gaia/star perpendicular control")
    field = {r["object_id"]: r["field"] for r in gaia}
    require(sum(r["mer_det_quality_flag"] == "0" and
                field[r["object_id"]] != "LDN1641" for r in selected) == 198,
            "clean quality-flag-zero subset")
    return len(selected), classes["weak or confused"]


def check_archive_and_rbh1():
    fields = rows("archive/tables/field_totals.csv")
    pairs = rows("archive/tables/archive_pairs.csv")
    hst = [r for r in fields if r["field"].startswith("xdf")]
    jwst = [r for r in fields if r["field"] == "JADES"]
    require(len(hst) == 5 and len(jwst) == 1, "archive field accounting")
    require(abs(sum(float(r["geometric_deg2"]) for r in hst) * 3600
                - 5.24288) < 1e-5, "blind HST area")
    require(abs(float(jwst[0]["geometric_deg2"]) * 3600
                - .26205) < 1e-5, "nested JWST area")
    require(tuple(total(hst, f) for f in
                  ("raw", "f814_raw_comparison", "pairs", "both_profile")) ==
            (637, 629, 104, 1), "HST dual-filter screen")
    require(len(pairs) == 107 and
            sum(r["field"] != "JADES" for r in pairs) == 104 and
            sum(r["both_profile"] == "True" for r in pairs) == 1,
            "matched archive pairs")
    require(tuple(int(jwst[0][f]) for f in
                  ("raw", "pairs", "both_profile")) == (38, 3, 0),
            "JWST screen")
    controls = blob("archive/tables/rbh1_measured.json")["ACS"]
    require({r["band"] for r in controls} == {"606", "814"} and
            all(r["within_one_arcsec"] is False for r in controls),
            "RBH-1: strict host misses in both bands")
    bands = {r["band"]: r for r in controls}
    require(round(bands["606"]["endpoint_to_host_arcsec"], 2) == 2.53 and
            round(bands["814"]["endpoint_to_host_arcsec"], 2) == 1.14 and
            bands["606"]["profile_gate"] is True and
            bands["814"]["profile_gate"] is False, "RBH-1 control metrics")


def main():
    a, b = check_coverage()
    n = check_recovery()
    c, u = check_catalog_and_reviews()
    check_archive_and_rbh1()
    print("PASS 700 receipts; 688 VIS/RMS coverage SHA-256 matches; frozen file digests")
    print(f"PASS 344 centres, 272 usable; footprints {a:.9f}, {b:.9f} deg2")
    print(f"PASS 8264/651 positive; 2529/450 negative; {n} injections")
    print(f"PASS 29953430 Q1 entries, {c} shape passes, {u} unresolved review rows")
    print("PASS HST/JWST and RBH-1 ledger")
    print("LIMIT: original FITS and full parent Parquet bytes were not rerun")


if __name__ == "__main__":
    main()
