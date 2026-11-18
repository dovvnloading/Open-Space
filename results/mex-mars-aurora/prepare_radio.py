"""Build reproducible one-second Mars Express Doppler pairs from native labels and bytes."""

import argparse
from collections import defaultdict
import csv
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
from statistics import median
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

from source import _columns, parse_rows

ROOT = Path(__file__).resolve().parent
RATIO = 3 / 11


def optional_median(values):
    valid = [value for value in values if value is not None]
    return median(valid) if valid else None


def reduce(rows):
    buckets = defaultdict(list)
    differences = []
    observed = [r["observed_hz"] for r in rows if r["observed_hz"] is not None and r["observed_hz"] > 0]
    for row in rows:
        obs, pred, residual = (row[key] for key in ("observed_hz", "predicted_hz", "residual_hz"))
        if obs is None or pred is None or residual is None or obs <= 0 or pred <= 0:
            continue
        differences.append(abs((obs - pred) - residual))
        buckets[row["utc"][:19]].append(row)
    per_second = {}
    for second, members in buckets.items():
        if len(members) < 5:
            continue
        per_second[second] = {
            "count": len(members),
            **{k: optional_median([row[k] for row in members]) for k in (
                "geometry_km", "observed_hz", "correction_hz", "residual_hz", "differential_hz", "sigma_hz"
            )},
        }
    return per_second, {
        "native_rows": len(rows),
        "positive_observed_rows": len(observed),
        "median_positive_native_observed_hz": median(observed) if observed else None,
        "missing_reported_residual_rows": sum(r["residual_hz"] is None for r in rows),
        "present_archive_differential_rows": sum(r["differential_hz"] is not None for r in rows),
        "valid_frequency_rows": len(differences),
        "eligible_seconds": len(per_second),
        "max_abs_residual_replay_error_hz": max(differences) if differences else None,
        "native_utc_first_last": [rows[0]["utc"], rows[-1]["utc"]],
    }


def groups_for(receipts):
    grouped = defaultdict(dict)
    for item in receipts:
        name = item["file"]
        if not name.endswith(".tab"):
            continue
        _, kind, _, _ = name.split("_")
        chain = name[3:7] if kind in ("d4s", "d5x") else "dsn"
        band = "s" if kind in ("d4s", "dps") else "x"
        group = grouped[(item["volume"], chain)]
        if band in group:
            raise ValueError(f"Repeated product for {group}")
        group[band] = name
    return grouped


def run(data_root, output, receipt_output):
    receipts = json.loads((ROOT / "radio-sources.json").read_bytes())
    for item in receipts:
        path = data_root / item["file"]
        if path.stat().st_size != item["bytes"] or sha256(path.read_bytes()).hexdigest() != item["raw_sha256"]:
            raise ValueError(f"Native archive product does not match receipt: {path}")
    rows = StringIO(newline="")
    fields = ("volume", "chain", "utc", "s_n", "x_n", "s_hz", "x_hz", "d_hz", "s_obs_hz", "x_obs_hz",
              "s_geometry_km", "x_geometry_km", "s_archive_differential_hz", "x_archive_differential_hz",
              "s_earth_correction_hz", "x_earth_correction_hz", "s_sigma_hz", "x_sigma_hz")
    writer = csv.writer(rows, lineterminator="\n")
    writer.writerow(fields)
    results = {}
    for (volume, chain), names in sorted(groups_for(receipts).items()):
        by_band = {}
        stats = {}
        for band in ("s", "x"):
            path = data_root / names[band]
            label = path.with_suffix(".lbl")
            by_band[band], stats[band] = reduce(parse_rows(path, label))
            stats[band]["observed_label_column"] = next(
                name for name in _columns(label.read_text(encoding="ascii")) if name.startswith("OBSERVED ")
            )
        s, x = by_band["s"], by_band["x"]
        common = sorted(s.keys() & x.keys())
        stats["common_seconds"] = len(common)
        if common:
            stats["common_start_stop"] = [common[0], common[-1]]
            stats["median_observed_hz"] = [median(band[t]["observed_hz"] for t in common) for band in (s, x)]
        for time in common:
            sl, xl = s[time], x[time]
            writer.writerow([
                volume, chain, time, sl["count"], xl["count"], sl["residual_hz"], xl["residual_hz"],
                sl["residual_hz"] - RATIO * xl["residual_hz"], sl["observed_hz"], xl["observed_hz"],
                sl["geometry_km"], xl["geometry_km"], sl["differential_hz"], xl["differential_hz"],
                sl["correction_hz"], xl["correction_hz"], sl["sigma_hz"], xl["sigma_hz"],
            ])
        results[f"{volume}-{chain}"] = stats
    with ZipFile(output, "w") as archive:
        info = ZipInfo("radio-seconds.csv")
        archive.writestr(info, rows.getvalue().encode("ascii"), compress_type=ZIP_DEFLATED, compresslevel=9)
    dsn = [stats for key, stats in results.items() if key.endswith("-dsn")]
    quality = {
        "released_dsn_passes": len(dsn),
        "dsn_s_files_labeled_x_band": sum(
            "X-BAND" in stats["s"]["observed_label_column"] for stats in dsn
        ),
        "dsn_s_files_carrier_2_to_3_ghz": sum(
            2e9 <= stats["s"]["median_positive_native_observed_hz"] <= 3e9 for stats in dsn
        ),
        "dsn_passes_with_calibrated_pairs": sum(stats["common_seconds"] > 0 for stats in dsn),
        "dsn_s_files_all_residuals_missing": sum(
            stats["s"]["missing_reported_residual_rows"] == stats["s"]["native_rows"] for stats in dsn
        ),
    }
    receipt_output.write_text(
        json.dumps({"quality": quality, "groups": results}, indent=2) + "\n",
        encoding="utf-8", newline="\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    run(args.data_root, args.output, args.receipt)
