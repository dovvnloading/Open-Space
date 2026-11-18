"""Measure paired agreement of the Sentinel-6 CoastWatch coastal products."""

import argparse
from datetime import datetime, timedelta, timezone
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from buoys import analyze_buoys

ROOT = Path(__file__).resolve().parent
EARTH_KM = 6371.0
EPOCH = datetime(1985, 1, 1, tzinfo=timezone.utc)
SECONDS_PER_DAY = 86400
FIRST_DAY = (datetime(2026, 9, 9, tzinfo=timezone.utc) - EPOCH).days
TRANSITION_DAY = (datetime(2026, 10, 30, tzinfo=timezone.utc) - EPOCH).days
HOLDOUT_DAY = (datetime(2026, 11, 1, tzinfo=timezone.utc) - EPOCH).days
END_DAY = (datetime(2026, 11, 16, tzinfo=timezone.utc) - EPOCH).days


def sphere(points):
    lat, lon = np.deg2rad(points).T
    return np.column_stack((np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)))


def arc(chord):
    return 2 * EARTH_KM * np.arcsin(np.minimum(chord / 2, 1))


def day_number(seconds):
    return np.floor_divide(seconds, SECONDS_PER_DAY).astype(int)


def day_label(day):
    return (EPOCH + timedelta(days=int(day))).date().isoformat()


def unique_pairs(a, b, max_km=5, symmetric=False, time_shift=0):
    a_xyz, b_xyz = sphere(a[:, 1:3]), sphere(b[:, 1:3])
    radius = 2 * np.sin(max_km / (2 * EARTH_KM))
    nearby = cKDTree(b_xyz).query_ball_point(a_xyz, radius)
    candidates = []
    for ai, options in enumerate(nearby):
        for bi in options:
            lag = b[bi, 0] - a[ai, 0] + time_shift
            effective_lag = abs(lag) if symmetric else lag
            if 10 <= effective_lag <= 70:
                separation = np.linalg.norm(a_xyz[ai] - b_xyz[bi])
                candidates.append((separation, abs(effective_lag - 30), ai, bi, lag))
    used_a, used_b, chosen = set(), set(), []
    for separation, _, ai, bi, lag in sorted(candidates):
        if ai in used_a or bi in used_b:
            continue
        used_a.add(ai)
        used_b.add(bi)
        chosen.append((ai, bi, arc(separation), lag))
    return np.array(sorted(chosen), dtype=float).reshape(-1, 4)


def numeric_summary(values):
    if not len(values):
        return {"n": 0, "median": None, "q25": None, "q75": None}
    q25, median, q75 = np.quantile(values, [0.25, 0.50, 0.75])
    return {"n": int(len(values)), "median": float(median), "q25": float(q25), "q75": float(q75)}


def day_bootstrap(visit_rows):
    days = sorted({row["date"] for row in visit_rows})
    if len(days) < 2:
        return None
    by_day = [np.array([row["near_minus_offshore_m"] for row in visit_rows if row["date"] == day])
              for day in days]
    sample = np.random.default_rng(6226).integers(len(days), size=(5000, len(days)))
    draws = [np.median(np.concatenate([by_day[index] for index in row])) for row in sample]
    return np.quantile(draws, [0.025, 0.975]).tolist()


def score(a, b, matches, coast):
    ai, bi = matches[:, :2].astype(int).T
    a, b = a[ai], b[bi]
    pair_time = (a[:, 0] + b[:, 0]) / 2
    coord = (a[:, 1:3] + b[:, 1:3]) / 2
    xyz = sphere(coord)
    gap = (np.diff(a[:, 0]) > 90) | (arc(np.linalg.norm(np.diff(xyz, axis=0), axis=1)) > 120)
    visit_id = np.cumsum(np.r_[False, gap])[:len(pair_time)]
    distance_to_coast = arc(cKDTree(sphere(coast)).query(xyz)[0])
    dates = day_number(pair_time)
    near = (distance_to_coast >= 20) & (distance_to_coast < 80)
    offshore = (distance_to_coast >= 150) & (distance_to_coast < 300)
    difference = b[:, 3] - a[:, 3]
    valid = np.isfinite(difference)
    wave = np.isfinite(a[:, 4]) & (a[:, 4] >= 0) & np.isfinite(b[:, 4]) & (b[:, 4] >= 0)
    average_wave = (a[:, 4] + b[:, 4]) / 2
    results = {"paired_coordinates": len(matches), "lag_s": numeric_summary(matches[:, 3]),
               "separation_km": numeric_summary(matches[:, 2]), "by_period": {}}
    visit_rows = []
    for visit in np.unique(visit_id):
        visit_rows_mask = visit_id == visit
        for day in np.unique(dates[visit_rows_mask]):
            group = visit_rows_mask & (dates == day)
            inshore = difference[group & near & valid]
            outshore = difference[group & offshore & valid]
            if len(inshore) >= 3 and len(outshore) >= 3:
                near_median = float(np.median(abs(inshore)))
                offshore_median = float(np.median(abs(outshore)))
                visit_rows.append({
                    "visit": int(visit), "date": day_label(day), "day": int(day),
                    "near_count": len(inshore), "offshore_count": len(outshore),
                    "near_median_abs_m": near_median, "offshore_median_abs_m": offshore_median,
                    "near_minus_offshore_m": near_median - offshore_median,
                    "near_signed_m": float(np.median(inshore)), "offshore_signed_m": float(np.median(outshore)),
                })
    for name, first, last in (
        ("reference", FIRST_DAY, TRANSITION_DAY), ("transition", TRANSITION_DAY, HOLDOUT_DAY),
        ("holdout", HOLDOUT_DAY, END_DAY),
    ):
        period = (dates >= first) & (dates < last)
        visits = [row for row in visit_rows if first <= row["day"] < last]
        results["by_period"][name] = {
            "coordinates": int(period.sum()), "valid_sla_pairs": int((period & valid).sum()),
            "valid_swh_pairs": int((period & wave).sum()),
            "coordinates_outside_primary_bins": int((period & ~near & ~offshore).sum()),
            "near_sla_difference_m": numeric_summary(difference[period & near & valid]),
            "offshore_sla_difference_m": numeric_summary(difference[period & offshore & valid]),
            "visit_near_minus_offshore_m": numeric_summary([item["near_minus_offshore_m"] for item in visits]),
            "visit_day_bootstrap_95_m": day_bootstrap(visits),
            "low_wave_abs_sla_difference_m": numeric_summary(abs(difference[period & valid & wave & (average_wave < 2)])),
            "high_wave_abs_sla_difference_m": numeric_summary(abs(difference[period & valid & wave & (average_wave >= 3)])),
        }
    results["supported_visits"] = visit_rows
    return results, (pair_time[wave], coord[wave], visit_id[wave], a[wave, 4], b[wave, 4])


def run(inputs, metadata):
    a, b, coast = (inputs[name] for name in ("6a", "6b", "coast"))
    main_pairs = unique_pairs(a, b)
    main, buoy_data = score(a, b, main_pairs, coast)
    tight, _ = score(a, b, unique_pairs(a, b, max_km=2), coast)
    symmetric, _ = score(a, b, unique_pairs(a, b, symmetric=True), coast)
    shifted = unique_pairs(a, b, time_shift=SECONDS_PER_DAY)
    mid = (a[main_pairs[:, 0].astype(int), 0] + b[main_pairs[:, 1].astype(int), 0]) / 2
    days = day_number(mid)
    daily = {}
    native = metadata["satellite"]
    for day in range(FIRST_DAY, END_DAY):
        label = day_label(day)
        stamp = label.replace("-", "")
        daily[label] = {
            instrument: native.get(f"{instrument}_{stamp}.nc", {}).get("native_rows")
            for instrument in ("6a", "6b")
        }
        daily[label]["regional_6a"] = int((day_number(a[:, 0]) == day).sum())
        daily[label]["regional_6b"] = int((day_number(b[:, 0]) == day).sum())
        selected = days == day
        daily[label]["coordinate_pairs"] = int(selected.sum())
        daily[label]["separation_km"] = numeric_summary(main_pairs[selected, 2])
        daily[label]["lag_s"] = numeric_summary(main_pairs[selected, 3])
    latency = {}
    for instrument in ("6a", "6b"):
        files = [row for filename, row in native.items() if filename.startswith(instrument)]
        latency[instrument] = {
            "files": len(files),
            "rows": sum(row["native_rows"] for row in files),
            "explicit_sla_latency_limit_files": sum(
                "-Qsla=latency" in row["history"] and "-Llatency=0,1" in row["history"] for row in files
            ),
        }
    return {
        "native_audit": latency,
        "daily_counts": daily,
        "primary": main,
        "two_km_sensitivity": tight,
        "symmetric_lag_sensitivity": symmetric,
        "one_day_shifted_coordinate_pairs": len(shifted),
        "buoys": analyze_buoys(inputs, *buoy_data),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    metadata = json.loads((ROOT / "native-metadata.json").read_bytes())
    path = ROOT / "regional-inputs.npz"
    if sha256(path.read_bytes()).hexdigest() != metadata["data_sha256"]:
        raise ValueError("Regional input digest mismatch")
    with np.load(path) as inputs:
        summary = run(inputs, metadata)
    args.output.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8", newline="\n")
