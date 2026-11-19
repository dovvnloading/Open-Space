#!/usr/bin/env python3
"""Fail-closed DJA v4.4 catalog-first z>9 selector; no FITS/astronomy results.

Requires exact publisher-MD5-verified 80,367-row v4.4 summary CSV.gz.
Does not download spectra. Each catalog row remains a separate spectrum.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

EXPECTED_ROWS = 80367
PUBLISHER_MD5 = "721c58c446f4a4eb1570aef02c4ab47b"
Z_COLS = ("z_best", "z_phot", "z", "zfit", "z_fit", "zsp_best",
          "z_16", "z_50", "z_84", "z16", "z50", "z84", "z_q84",
          "z_p84", "z_upper", "z_hi", "z_median", "z_best_hi")
HERE = Path(__file__).resolve().parent
CONTROLS = HERE.parent / "jwst-dja-nirspec-highz-screen" / "literature_baseline.json"


def hashes(path):
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            sha.update(chunk); md5.update(chunk)
    return sha.hexdigest(), md5.hexdigest()


def angular_sep_arcsec(a, b, c, d):
    r1, d1, r2, d2 = map(math.radians, (a, b, c, d))
    hav = math.sin((d2-d1)/2)**2 + math.cos(d1)*math.cos(d2)*math.sin((r2-r1)/2)**2
    return 2 * math.asin(min(1., math.sqrt(max(0., hav)))) * 180/math.pi * 3600


def classify(row, zcols, controls, radius=0.5):
    hints = []
    for col in zcols:
        try:
            val = float(row.get(col, ""))
            if math.isfinite(val) and val > 9:
                hints.append(f"{col}={val:.8g}")
        except (ValueError, TypeError):
            pass
    matches = []
    file = row.get("file", "")
    for ctrl in controls:
        key = ctrl["catalog_source_id"]
        id_match = file.endswith("_" + key + ".spec.fits")
        try:
            sky_match = angular_sep_arcsec(float(row["ra"]), float(row["dec"]),
                                          float(ctrl["ra_deg"]), float(ctrl["dec_deg"])) <= radius
        except (ValueError, TypeError, KeyError):
            sky_match = False
        if id_match or sky_match:
            matches.append(f"{ctrl['name']}:{'id' if id_match else 'position'}")
    return hints, matches


def select(catalog, out, controls_path=CONTROLS, expected=EXPECTED_ROWS,
           expected_md5=PUBLISHER_MD5):
    out.mkdir(parents=True, exist_ok=True)
    sha, md5 = hashes(catalog)
    if md5.lower() != expected_md5.lower():
        raise ValueError(f"DJA v4.4 publisher MD5 mismatch: observed={md5}, expected={expected_md5}")
    controls = json.loads(controls_path.read_text())["published_positive_controls"]
    tmp_audit = out / "catalog_audit.csv.partial"
    tmp_selected = out / "selected_spectra.csv.partial"
    total = selected = 0; matched = set(); zcols = []
    try:
        with gzip.open(catalog, "rt", encoding="utf-8", newline="") as inp, \
             tmp_audit.open("w", encoding="utf-8", newline="") as audit, \
             tmp_selected.open("w", encoding="utf-8", newline="") as targets:
            reader = csv.DictReader(inp)
            cols = reader.fieldnames or []
            for col in ("root", "file", "ra", "dec", "srcid", "z_best", "z_phot"):
                if col not in cols:
                    raise ValueError(f"Missing DJA catalog field: {col}")
            zcols = [z for z in Z_COLS if z in cols]
            fields = cols + ["dja_catalog_row_index", "targeted_selected", "highz_hints", "published_controls"]
            a = csv.DictWriter(audit, fieldnames=fields); s = csv.DictWriter(targets, fieldnames=fields)
            a.writeheader(); s.writeheader()
            for index, row in enumerate(reader):
                hints, matches = classify(row, zcols, controls)
                keep = bool(hints or matches)
                row.update({"dja_catalog_row_index": index, "targeted_selected": str(keep).lower(),
                            "highz_hints": ";".join(hints), "published_controls": ";".join(matches)})
                a.writerow(row)
                if keep:
                    root, file = row["root"], row["file"]
                    if any(p in ("", ".", "..") for p in (root, file)) or "/" in root + file or "\\" in root + file or not file.endswith(".spec.fits"):
                        raise ValueError(f"Unsafe/unexpected DJA spectrum identifier at row {index}")
                    s.writerow(row); selected += 1
                    matched.update(m.split(":")[0] for m in matches)
                total += 1
        if total != expected:
            raise ValueError(f"Incomplete DJA v4.4 summary: {total} rows, expected {expected}")
        audit_path, selection_path = out / "catalog_audit.csv", out / "selected_spectra.csv"
        tmp_audit.replace(audit_path); tmp_selected.replace(selection_path)
        receipt = {"scope": "DJA v4.4 catalog hints only; NO independent quick fits",
                   "catalog_sha256": sha, "catalog_md5": md5, "catalog_rows": total,
                   "selected_spectrum_rows": selected, "z_hint_fields": zcols,
                   "published_controls_matched": sorted(matched),
                   "published_controls_unmatched": sorted(set(c["name"] for c in controls)-matched),
                   "selection_sha256": hashes(selection_path)[0], "audit_sha256": hashes(audit_path)[0],
                   "status": "SELECTION ONLY; no FITS downloads, 1D fits, 2D checks or NIRCam photometry",
                   "selection_conservatism": "All rows with any z hint >9 retained, irrespective of grade; controls by exact filename suffix or <=0.5 arcsec. Untested z>9 solutions outside catalog hints cannot be excluded."}
        (out / "selection_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
        return receipt
    finally:
        tmp_audit.unlink(missing_ok=True); tmp_selected.unlink(missing_ok=True)


def self_test():
    controls = [{"name": "control", "catalog_source_id": "999_5", "ra_deg": 10., "dec_deg": 20.}]
    rows = [
        {"z_best": "9.2", "z_phot": "", "file": "r_prism-clear_1_1.spec.fits", "ra": "0", "dec": "0"},
        {"z_best": "4", "z_phot": "9.01", "file": "r_prism-clear_1_2.spec.fits", "ra": "0", "dec": "0"},
        {"z_best": "5", "z_phot": "NaN", "file": "r_prism-clear_1_3.spec.fits", "ra": "10", "dec": "20"},
        {"z_best": "5", "z_phot": "7", "file": "r_prism-clear_1_4.spec.fits", "ra": "0", "dec": "0"},
        {"z_best": "5", "z_phot": "7", "file": "r_prism-clear_999_5.spec.fits", "ra": "0", "dec": "0"},
    ]
    reasons = [classify(r, ("z_best", "z_phot"), controls) for r in rows]
    assert [bool(a or b) for a,b in reasons] == [True, True, True, False, True]
    assert reasons[2][1] == ["control:position"]
    assert reasons[4][1] == ["control:id"]
    assert angular_sep_arcsec(10,20,10,20)==0
    print("PASS selection rule / positive-control positional & ID inclusion (synthetic only)")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--catalog", type=Path)
    p.add_argument("--out", type=Path)
    p.add_argument("--self-test", action="store_true")
    args=p.parse_args()
    if args.self_test: self_test()
    else:
        if not args.catalog or not args.out: p.error("--catalog and --out required")
        print(json.dumps(select(args.catalog, args.out), indent=2))
