"""Independent archival checks and expanded CORALIE epochs for HD 25673."""

import argparse
from collections import Counter
import csv
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
from scipy import stats

from screen import HERE, clean_number, finite, linear, yearly


def records(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        if stream.read(1) == "#":
            stream.readline()
        else:
            stream.seek(0)
        return list(csv.DictReader(stream))


def wls(x, y, error):
    x = np.asarray(x) / 365.25
    weights = 1.0 / np.asarray(error) ** 2
    x -= np.average(x, weights=weights)
    centered_y = y - np.average(y, weights=weights)
    slope_variance = 1 / np.dot(weights, x ** 2)
    slope = float(slope_variance * np.dot(weights * x, centered_y))
    residual = centered_y - slope * x
    chisq = float(np.dot(weights, residual ** 2))
    return {"slope_mps_per_year": slope, "formal_error_mps_per_year": float(np.sqrt(slope_variance)),
            "chi2": chisq, "dof": len(x) - 2,
            "weighted_rms_mps": float(np.sqrt(np.average(residual ** 2, weights=weights)))}


def describe_coralie(rows, fp_only):
    accepted = [item for item in rows if item["outcome"] == "ok" and finite(item, "qc") == 1
                and item["mask"] and 0 < finite(item, "errmps") <= 20 and
                np.isfinite(finite(item, "rvmps")) and
                item["object"].replace(" ", "").upper() == "HD25673"
                and (not fp_only or item["ref"] == "OBJ_FP")]
    per_mask = Counter(item["mask"] for item in accepted)
    if not per_mask:
        return {"usable": 0}
    chosen = min(per_mask, key=lambda name: (-per_mask[name], name))
    accepted = [item for item in accepted if item["mask"] == chosen]
    annual = yearly(accepted, "rvmps")
    slope, p, loo = linear(annual, 0) if len(annual) >= 3 else (None, None, [])
    return {"usable": len(accepted), "mask": chosen,
            "reference": dict(Counter(row["ref"] for row in accepted)), "yearly": annual,
            "slope_mps_per_year": slope, "two_sided_p": p, "leave_one_year_out": list(map(float, loo))}


def replay():
    harps_file = HERE / "harps-rvbank.csv"
    expected = "90107e3d0ffc66e50033d5a53e74620d102993819359301f0f08f262f1e3daa8"
    with harps_file.open("rb") as stream:
        if file_digest(stream, "sha256").hexdigest() != expected:
            raise ValueError("HARPS RVBank bytes differ from the saved response")
    harps = [row for row in records(harps_file) if int(row["f_RV"]) == 0]
    obs = [dict(row, mjd=float(row["BJD"]) - 2400000.5, year=int(row["TIMEID"][:4])) for row in harps]
    epochs = np.array([float(row["BJD"]) - float(harps[0]["BJD"]) for row in harps])
    harps_result = {}
    for rv, error in (("RV_mlc_nzp", "e_RV_mlc_nzp"), ("RV_drs_nzp", "e_RV_drs_nzp")):
        v = np.array([float(item[rv]) for item in harps])
        result = wls(epochs, v, [float(item[error]) for item in harps])
        bins = yearly(obs, rv)
        slope, p, loo = linear(bins, 0)
        yearfit = {"slope_mps_per_year": slope, "two_sided_p": p, "leave_one_year_out": list(map(float, loo))}
        residual = v - result["slope_mps_per_year"] * epochs / 365.25
        ranks = {}
        for key in ("Halpha", "RHKp", "BIS"):
            values = np.array([finite(row, key) for row in harps])
            ok = np.isfinite(values) & (values > -999)
            ranks[key] = clean_number(stats.spearmanr(residual[ok], values[ok]).statistic) if ok.sum() >= 8 else None
        harps_result[rv] = {"weighted_fit": result, "yearly": bins, "yearly_fit": yearfit,
                            "detrended_activity_rank_correlation": ranks}
    coralie = records(HERE / "hd25673-headers.csv")
    phases = ("coralie98", "coralie07", "training", "holdout")
    expanded = {}
    for phase in phases:
        subset = [row for row in coralie if row["split"] == phase]
        entry = describe_coralie(subset, False)
        entry["selected"] = len(subset)
        if phase in ("training", "holdout"):
            entry["fabry_perot_only"] = describe_coralie(subset, True)
        expanded[phase] = entry

    earlier_coralie = expanded["coralie07"]["yearly"]
    earlier_harps = harps_result["RV_mlc_nzp"]["yearly"]
    common_years = sorted({row["year"] for row in earlier_coralie} & {row["year"] for row in earlier_harps})
    overlap = {
        "years": common_years,
        "coralie07_slope_mps_per_year": linear([row for row in earlier_coralie if row["year"] in common_years], 0)[0],
        "harps_slope_mps_per_year": linear([row for row in earlier_harps if row["year"] in common_years], 0)[0],
    }
    modern = expanded["training"]["fabry_perot_only"]["yearly"] + expanded["holdout"]["fabry_perot_only"]["yearly"]
    modern_slope = linear(modern, 0)[0]
    apogee = records(HERE / "apogee-visits.csv")
    apo_epochs = [float(row["jd"]) for row in apogee]
    return {
        "coralie_all_nights": expanded,
        "harps_coralie07_same_years": overlap,
        "combined_modern_fp": {"years": len(modern), "slope_mps_per_year": modern_slope,
                               "first_to_last_year_median_change_mps": modern[-1]["value"] - modern[0]["value"]},
        "harps": {"observations": len(harps), "span_days": float(np.ptp(epochs)), "fits": harps_result},
        "apogee": {"observations": len(apogee), "span_days": float(np.ptp(apo_epochs)),
                   "rms_heliocentric_rv_mps": float(1000 * np.std([float(r["vhelio"]) for r in apogee]))},
    }


def plot(result, path):
    import matplotlib.pyplot as plt

    entries = [
        ("CORALIE98", "#896a40", result["coralie_all_nights"]["coralie98"]["yearly"]),
        ("CORALIE07", "#8559a5", result["coralie_all_nights"]["coralie07"]["yearly"]),
        ("HARPS", "#3c527d", result["harps"]["fits"]["RV_mlc_nzp"]["yearly"]),
        ("CORALIE14 training", "#127a75", result["coralie_all_nights"]["training"]["yearly"]),
        ("CORALIE14 holdout", "#a13d49", result["coralie_all_nights"]["holdout"]["yearly"]),
    ]
    fig, axes = plt.subplots(2, 1, figsize=(9, 6.1), sharex=True, layout="constrained")
    for index, (name, color, annual) in enumerate(entries):
        ax = axes[int(index > 2)]
        x = [2000 + (row["mjd"] - 51544.5) / 365.25 for row in annual]
        y = np.array([row["value"] for row in annual])
        ax.plot(x, y - np.median(y), "o-", color=color, label=name, linewidth=1.6)
    axes[0].set_title("HD 25673’s line-of-sight velocity fell in both archives", loc="left", fontsize=14)
    for ax in axes:
        ax.set_ylabel("Velocity relative to\nseries median (m/s)")
        ax.axhline(0, color="#c2c7cd", linewidth=1)
        ax.grid(axis="y", color="#e5e8eb")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(loc="lower left", frameon=False, ncols=2, fontsize=8)
    axes[1].set_xlabel("Observation year")
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()
    result = replay()
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    if args.plot:
        plot(result, args.plot)


if __name__ == "__main__":
    main()
