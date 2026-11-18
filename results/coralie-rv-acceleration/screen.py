"""Replay the frozen two-epoch velocity screen and identify HGCA counterparts."""

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import file_digest
from itertools import islice, permutations
import json
import math
from pathlib import Path
import re

from astropy.io import fits
import numpy as np
from scipy import stats


HERE = Path(__file__).resolve().parent
YEAR_DAYS = 365.25
SPLITS = ("training", "holdout")
CONTROLS = {
    "HD10700": "quiet",
    "HD20794": "quiet",
    "HD219077": "known companion",
    "HD211847": "known companion",
}


def finite(row, key):
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError):
        return math.nan
    return value if math.isfinite(value) else math.nan


def clean_number(value, digits=6):
    return round(float(value), digits) if math.isfinite(value) else None


def yearly(data, quantity):
    bins = defaultdict(list)
    for item in data:
        if math.isfinite(finite(item, quantity)):
            bins[int(item["year"])].append(item)
    return [{"year": year, "mjd": float(np.median([float(r["mjd"]) for r in rows])),
             "value": float(np.median([float(r[quantity]) for r in rows])), "nights": len(rows)}
            for year, rows in sorted(bins.items())]


def slope_weights(dates):
    x = (np.asarray(dates) - np.mean(dates)) / YEAR_DAYS
    return x / np.dot(x, x)


def linear(bins, tail):
    dates = np.array([v["mjd"] for v in bins])
    actual = np.array([v["value"] for v in bins])
    weights = slope_weights(dates)
    slope = np.dot(weights, actual)
    rearrangements = permutations(actual)
    exceedances = 0
    for chunk in iter(lambda: list(islice(rearrangements, 8192)), []):
        permuted = np.array(chunk) @ weights
        if tail == 0:
            exceedances += np.count_nonzero(abs(permuted) >= abs(slope) - 1e-9)
        else:
            exceedances += np.count_nonzero(tail * permuted >= tail * slope - 1e-9)
    p = exceedances / math.factorial(len(actual))
    loo = [np.dot(slope_weights(np.delete(dates, i)), np.delete(actual, i)) for i in range(len(bins))]
    return float(slope), float(p), loo


def indicator(rows, key, slope):
    bins = yearly(rows, key)
    shape_slope = linear(bins, 0)[0] if len(bins) >= 4 else math.nan
    jointly_measured = [row for row in rows if math.isfinite(finite(row, key))]
    rv = np.array([finite(row, "rvmps") for row in jointly_measured])
    x = np.array([float(row["mjd"]) for row in jointly_measured])
    shape = [finite(row, key) for row in jointly_measured]
    residual = rv - slope * (x - x.mean()) / YEAR_DAYS if len(x) else rv
    rank = stats.spearmanr(residual, shape).statistic if (
        len(rv) >= 8 and np.ptp(residual) > 0 and np.ptp(shape) > 0
    ) else math.nan
    return {"nights": len(rv), "years": len(bins), "shape_slope_mps_per_year": clean_number(shape_slope),
            "detrended_rv_rank_correlation": clean_number(rank)}


def fit_target(rows):
    eligible = []
    quality_failures = Counter()
    for row in rows:
        reason = None
        if row["outcome"] != "ok":
            reason = "header_error"
        elif finite(row, "qc") != 1:
            reason = "science_qc"
        elif not (math.isfinite(finite(row, "rvmps")) and 0 < finite(row, "errmps") <= 20 and row["mask"]):
            reason = "velocity_error_or_mask"
        elif re.sub(r"\s+", "", row["target"].strip().upper()) != re.sub(r"\s+", "", row["object"].strip().upper()):
            reason = "target_label"
        if reason:
            quality_failures[reason] += 1
        else:
            eligible.append(row)
    masks = Counter(row["mask"] for row in eligible)
    mask = min(masks, key=lambda name: (-masks[name], name), default=None)
    phases = {s: [row for row in eligible if row["split"] == s and row["mask"] == mask] for s in SPLITS}
    bins = {s: yearly(value, "rvmps") for s, value in phases.items()}
    summary = {
        "selected_nights": {s: sum(row["split"] == s for row in rows) for s in SPLITS},
        "mask": mask, "acceptable_masks": dict(masks), "quality_exclusions": dict(quality_failures),
        "usable_nights": {s: len(r) for s, r in phases.items()},
        "years": {s: len(b) for s, b in bins.items()},
        "span_days": {s: round(max(float(r["mjd"]) for r in rs) - min(float(r["mjd"]) for r in rs), 3)
                      if rs else 0 for s, rs in phases.items()},
        "calibration": {s: dict(Counter(f"{r['ref']}|{r['fp']}" for r in rows)) for s, rows in phases.items()},
        "yearly": bins,
    }
    summary["fit_available"] = all(summary["usable_nights"][s] >= 10 and summary["years"][s] >= 4
                                   and summary["span_days"][s] >= 700 for s in SPLITS)
    summary["passes"] = False
    if not summary["fit_available"]:
        return summary
    first = linear(bins["training"], 0)
    sign = int(np.sign(first[0]))
    second = linear(bins["holdout"], sign)
    diagnostics = {}
    for phase, (slope, p, loo) in zip(SPLITS, (first, second)):
        diagnostics[phase] = {
            "slope_mps_per_year": clean_number(slope), "permutation_p": p,
            "loo_slopes_mps_per_year": list(map(clean_number, loo)),
            "bisector": indicator(phases[phase], "bismps", slope),
            "fwhm": indicator(phases[phase], "fwhmmps", slope),
        }
    summary["fits"] = diagnostics
    summary["passes"] = all(sign * slope >= 10 and p <= 0.05 and all(sign * x > 0 for x in loo)
                            for slope, p, loo in (first, second))
    summary["fisher_p"] = float(stats.chi2.sf(-2 * math.log(first[1] * second[1]), 4))
    return summary


def controls():
    rows = []
    for phase, start in (("training", 2015), ("holdout", 2020)):
        for i in range(5):
            for offset in (0, 30):
                rows.append({
                    "target": "TEST", "object": "TEST", "split": phase, "year": start + i,
                    "mjd": 57260 + (start - 2015 + i) * YEAR_DAYS + offset,
                    "rvmps": 20 * (i + offset / YEAR_DAYS), "errmps": 2,
                    "outcome": "ok", "mask": "K5", "qc": 1, "ref": "FP", "fp": "FP",
                })
    consistent = fit_target(rows)["passes"]
    reversal = [dict(r, rvmps=-r["rvmps"]) if r["split"] == "holdout" else r for r in rows]
    reversed_pass = fit_target(reversal)["passes"]
    if not consistent or reversed_pass:
        raise RuntimeError("Synthetic prediction controls failed")
    return {"positive_linear_passed": consistent, "sign_reversal_passed": reversed_pass}


def typical_point(rows):
    radians = np.deg2rad([float(row["ra"]) for row in rows])
    origin = np.arctan2(np.sin(radians).mean(), np.cos(radians).mean())
    ra = np.rad2deg(origin + np.median(np.angle(np.exp(1j * (radians - origin))))) % 360
    dec = np.median([float(row["dec"]) for row in rows])
    return float(ra), float(dec)


def sky_offsets(catalog, ra, dec):
    dr = np.deg2rad(catalog["gaia_ra"] - ra)
    dd = np.deg2rad(catalog["gaia_dec"] - dec)
    half = np.sin(dd / 2) ** 2 + np.cos(np.deg2rad(dec)) * np.cos(np.deg2rad(catalog["gaia_dec"])) * np.sin(dr / 2) ** 2
    return 7200 * np.rad2deg(np.arcsin(np.sqrt(np.clip(half, 0, 1))))


def asymmetric_pm(catalog, row, separation):
    point = catalog[row]
    chi2 = float(point["chisq"])
    pm_norm = np.hypot(point["pmra_gaia"], point["pmdec_gaia"])
    # RV change from the rotation of the line of sight, using distance from parallax.
    pc_m = 3.085677581491367e16
    seconds_per_year = YEAR_DAYS * 86400
    mu_per_year = np.deg2rad(pm_norm / 3.6e6)
    parallax = point["parallax_gaia"]
    perspective = mu_per_year ** 2 * 1000 * pc_m / (parallax * seconds_per_year) if parallax > 0 else math.nan
    return {
        "hip": int(point["hip_id"]), "gaia_edr3": str(point["gaia_source_id"]),
        "separation_arcsec": clean_number(separation), "gaia_ra_deg": clean_number(point["gaia_ra"], 8),
        "gaia_dec_deg": clean_number(point["gaia_dec"], 8), "constant_pm_chi2_2dof": clean_number(chi2),
        "constant_pm_nominal_p": float(stats.chi2.sf(chi2, 2)),
        "parallax_mas": clean_number(parallax), "gaia_pm_mas_per_year": clean_number(pm_norm),
        "perspective_rv_acceleration_mps_per_year": clean_number(perspective),
        "hgca_rv_kms": clean_number(point["radial_velocity"]),
        "hgca_rv_source": str(point["radial_velocity_source"]).strip(),
        "nonlinear_dpmra": clean_number(point["nonlinear_dpmra"]),
        "nonlinear_dpmdec": clean_number(point["nonlinear_dpmdec"]),
    }


def match_astrometry(catalog, rows):
    ra, dec = typical_point(rows)
    separations = sky_offsets(catalog, ra, dec)
    nearest = np.flatnonzero(separations <= 30)
    identifier = re.fullmatch(r"HIP\s*(\d+)", rows[0]["target"])
    match = "no match"
    candidates = nearest
    compact_id = None
    if identifier:
        hip = int(identifier.group(1))
        hits = np.flatnonzero(catalog["hip_id"] == hip)
        candidates = np.union1d(nearest, hits)
        if len(hits) == 1 and separations[hits[0]] <= 30:
            compact_id = hits[0]
            match = "HIP and position"
        else:
            match = "unverified HIP identifier"
    elif len(nearest) == 1:
        compact_id = nearest[0]
        match = "single positional"
    elif len(nearest) > 1:
        match = "ambiguous position"
    return {"target_median_radec": [clean_number(ra, 8), clean_number(dec, 8)], "status": match,
            "within_30_arcsec": len(nearest),
            "candidates": [asymmetric_pm(catalog, i, separations[i]) for i in candidates],
            "selected": None if compact_id is None else asymmetric_pm(catalog, compact_id, separations[compact_id])}


def analyze(observations, astrometry):
    synthetic = controls()
    receipt = json.loads((HERE / "hgca-receipt.json").read_text(encoding="utf-8"))
    with astrometry.open("rb") as stream:
        digest = file_digest(stream, "sha256").hexdigest()
    if digest != receipt["sha256"]:
        raise ValueError("HGCA FITS checksum mismatch")
    grouped = defaultdict(list)
    with observations.open(encoding="utf-8", newline="") as stream:
        for record in csv.DictReader(stream):
            grouped[record["target"]].append(record)
    with fits.open(astrometry) as product:
        catalog = product[1].data
        all_targets = []
        for target, rows in sorted(grouped.items()):
            fitted = fit_target(rows)
            astro = match_astrometry(catalog, rows) if fitted["passes"] else None
            all_targets.append(dict(fitted, target=target, astrometry=astro, benchmark=CONTROLS.get(target.replace(" ", ""))))
    tested = sorted((row for row in all_targets if row["fit_available"]), key=lambda row: row["fisher_p"])
    n = len(tested)
    smallest = 1
    for index in range(n - 1, -1, -1):
        row = tested[index]
        smallest = min(smallest, n * row["fisher_p"] / (index + 1))
        row["bh_q"] = round(smallest, 8)
    phases = {phase: [row["fits"][phase]["slope_mps_per_year"] for row in tested] for phase in SPLITS}
    passing = [row for row in tested if row["passes"]]
    notable = [row for row in passing if row["bh_q"] <= 0.1]
    return {
        "summary": {
            "selected_spectra": sum(map(len, grouped.values())),
            "available_headers": sum(row["outcome"] == "ok" for group in grouped.values() for row in group),
            "selected_targets": len(grouped), "tested_targets": n,
            "passed_two_epoch_gate": len(passing), "passed_two_epoch_and_bh": len(notable),
            "sample_median_slopes_mps_per_year": {k: clean_number(np.median(v)) if v else None for k, v in phases.items()},
            "passing_targets": [row["target"] for row in passing],
        },
        "controls": synthetic, "targets": all_targets,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("observations", type=Path)
    parser.add_argument("astrometry", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    results = analyze(args.observations, args.astrometry)
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(results["summary"], indent=2))


if __name__ == "__main__":
    main()
