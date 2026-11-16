#!/usr/bin/env python3
"""External-institution stress test: Goddard GMM-3 versus two JPL Mars models.

GMM-3 uses its own NOT-IAU body-fixed frame and partly overlapping tracking arcs.
Results are *not* independent evidence for temporal variation in Mars gravity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import re

import analyze as g

GMM_TABLE = "gmm3_120_sha.tab"
GMM_LABEL = "gmm3_120_sha.lbl"


def read_gmm3(cache: pathlib.Path) -> tuple[dict, dict]:
    binary = g.download(g.ARCHIVE + GMM_TABLE, cache / GMM_TABLE)
    label = g.download(g.ARCHIVE + GMM_LABEL, cache / GMM_LABEL)
    txt = label.decode("ascii").upper()
    assert re.search(r'ORIGINAL_PRODUCT_ID\s*=\s*"GMM-3"', txt)
    assert re.search(r"ROWS\s*=\s*7380\b", txt)
    assert "NOT* THE IAU" in txt or "NOT THE IAU" in txt
    parsed = g.parse_table(binary)
    report = {
        "table_url": g.ARCHIVE + GMM_TABLE,
        "label_url": g.ARCHIVE + GMM_LABEL,
        "sha256_table": hashlib.sha256(binary).hexdigest(),
        "sha256_label": hashlib.sha256(label).hexdigest(),
        "table_bytes": len(binary),
        "label_bytes": len(label),
        "rows": len(parsed["coefficients"]),
        "GM_km3_s2": parsed["GM_km3_s2"],
        "reference_warning": "GMM-3 PDS label: reference system is NOT the IAU Mars system",
    }
    return parsed, report


def run(cache: pathlib.Path, out: pathlib.Path) -> dict:
    models = {}
    sources = {}
    for model, name in g.FILES.items():
        raw = g.download(g.ARCHIVE + name, cache / name)
        label = g.download(g.ARCHIVE + g.LABELS[model], cache / g.LABELS[model])
        assert hashlib.sha256(raw).hexdigest() == g.EXPECTED_SHA256[model]["table"]
        assert hashlib.sha256(label).hexdigest() == g.EXPECTED_SHA256[model]["label"]
        g.validate_label(label, model)
        models[model] = g.parse_table(raw)["coefficients"]
        sources[model] = {"sha256_table": hashlib.sha256(raw).hexdigest()}
    independent, rec = read_gmm3(cache)
    models["GMM3"] = independent["coefficients"]
    sources["GMM3"] = rec
    results = {}
    for source, target in [("GMM3", "MRO120F"),
                           ("GMM3", "MRO120G"),
                           ("MRO120F", "MRO120G")]:
        a, b = models[source], models[target]
        yaw = g.fit_yaw(a, b, 2, 10)
        results[source + "_to_" + target] = {
            "fit_training_degrees": [2, 10],
            "effective_yaw_arcsec": yaw * g.ARCSEC_PER_RAD,
            "heldout_individually_fitted_yaw_arcsec": {
                "11-20": g.fit_yaw(a, b, 11, 20) * g.ARCSEC_PER_RAD,
                "21-40": g.fit_yaw(a, b, 21, 40) * g.ARCSEC_PER_RAD,
            },
            "trained_fraction_mismatch_removed": g.band_stats(
                a, b, 2, 10, yaw, 42828.375, 400
            )["fraction_raw_power_explained_by_rotation"],
            "heldout_11_20_fraction_removed": g.band_stats(
                a, b, 11, 20, yaw, 42828.375, 400
            )["fraction_raw_power_explained_by_rotation"],
            "heldout_21_40_fraction_removed": g.band_stats(
                a, b, 21, 40, yaw, 42828.375, 400
            )["fraction_raw_power_explained_by_rotation"],
        }
    triangle = (results["GMM3_to_MRO120G"]["effective_yaw_arcsec"]
                - results["GMM3_to_MRO120F"]["effective_yaw_arcsec"]
                - results["MRO120F_to_MRO120G"]["effective_yaw_arcsec"])
    report = {
        "purpose": "Stress test of cross-institution frame sensitivity; not a fit to new tracking data",
        "sources": sources, "pairs": results,
        "pairwise_angle_triangle_nonclosure_arcsec": triangle,
        "interpretation_warning": [
            "GSFC GMM3 (2016) is independently inverted but uses some same MRO/MGS/Odyssey tracking as JPL models.",
            "GMM3's label explicitly says its reference frame is NOT IAU; its fitted yaw has no direct absolute geographic meaning.",
            "Pairwise best-fit yaw is not an exactly transitive physical rotation when model harmonics change.",
            "Differences can reflect model data epochs, different regularization, tides, reference/pole, or data weighting.",
            "Do not report these differences as secular variation in Mars mass, mantle or rotation."
        ]
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"pairs": results, "triangle_nonclosure_arcsec": triangle}, indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=pathlib.Path, default=pathlib.Path("/tmp/open-space-mro120g-pds"))
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("/tmp/mars-gravity-gmm3-control.json"))
    opts = parser.parse_args()
    run(opts.cache, opts.output)
