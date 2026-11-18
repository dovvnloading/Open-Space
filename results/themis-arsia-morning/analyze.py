"""Replay the frozen Arsia matched-footprint brightening screen."""

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SETS = {"search": (-10.0, -8.5), "south": (-10.55, -10.05), "north": (-8.45, -7.95)}
BLOCK_DEG = 0.05


def interval(values, name):
    low, high = SETS[name]
    return (values >= low) & (values < high)


def available(pair):
    return np.isfinite(pair["earlier"]).all(axis=0) & np.isfinite(pair["later"]).all(axis=0)


def longest_run(flags):
    positions = np.flatnonzero(flags)
    lengths = [len(group) for group in np.split(positions, np.flatnonzero(np.diff(positions) != 1) + 1)]
    return max(lengths, default=0)


def profile(pair, constant=False):
    lat = np.broadcast_to(pair["lat_deg"][:, None], pair["earlier"].shape[1:])
    shared = available(pair)
    latitude = lat[shared]
    log_ratio = (np.log(pair["later"][:, shared]) - np.log(pair["earlier"][:, shared])).T
    shoulders = [interval(latitude, "south"), interval(latitude, "north")]
    counts = {name: int(interval(latitude, name).sum()) for name in SETS}
    result = {
        "paired_cells": int(shared.sum()), "cells_in_intervals": counts,
        "raw_search_median_change": None, "corrected_search_median_change": None,
        "shoulder_peak_positive_change": None, "longest_coincident_run_blocks": None,
        "passes_gate": False, "blocks": [],
    }
    if min(counts["south"], counts["north"]) < 30:
        return result
    if min(counts["search"], len(log_ratio)) == 0:
        return result
    if constant:
        background = np.median(log_ratio[shoulders[0] | shoulders[1]], axis=0)
    else:
        lower, upper = [np.median(log_ratio[mask], axis=0) for mask in shoulders]
        south, north = [np.median(latitude[mask]) for mask in shoulders]
        background = lower + ((latitude - south) / (north - south))[:, None] * (upper - lower)
    change = np.expm1(log_ratio - background)
    assigned = np.floor(latitude / BLOCK_DEG).astype(int)
    blocks = np.arange(np.floor(-10.65 / BLOCK_DEG), np.ceil(-7.85 / BLOCK_DEG), dtype=int)
    midpoint = BLOCK_DEG * (blocks + 0.5)
    medians = np.full((len(blocks), log_ratio.shape[1]), np.nan)
    block_counts = []
    for i, index in enumerate(blocks):
        selected = assigned == index
        n = int(selected.sum())
        block_counts.append(n)
        if n >= 10:
            medians[i] = np.median(change[selected], axis=0)
    valid_blocks = np.isfinite(medians).all(axis=1)
    control = valid_blocks & (interval(midpoint, "south") | interval(midpoint, "north"))
    search = interval(latitude, "search")
    result["raw_search_median_change"] = np.expm1(np.median(log_ratio[search], axis=0)).tolist()
    result["corrected_search_median_change"] = np.median(change[search], axis=0).tolist()
    result["blocks"] = [
        {"latitude_deg": float(lat_mid), "common_cells": n,
         "median_change": row.tolist() if valid else None}
        for lat_mid, n, row, valid in zip(midpoint, block_counts, medians, valid_blocks)
    ]
    if not control.any():
        return result
    peak = np.maximum(0, medians[control].max(axis=0))
    signal = interval(midpoint, "search") & valid_blocks
    signal &= (medians > np.maximum(0.10, peak)).all(axis=1)
    run = longest_run(signal)
    result["shoulder_peak_positive_change"] = peak.tolist()
    result["longest_coincident_run_blocks"] = run
    result["passes_gate"] = run >= 6
    return result


def load_pair(bundle, name, scale):
    return {key: bundle[f"{name}_{scale}_{key}"] for key in
            ("east_m", "north_m", "lat_deg", "earlier", "later")}


def shift(pair, dimension, delta):
    later = np.roll(pair["later"], delta, axis=dimension)
    edge = 0 if delta > 0 else -1
    if dimension == 2:
        later[:, :, edge] = np.nan
    else:
        later[:, edge, :] = np.nan
    return {**pair, "later": later}


def fake_patch(pair):
    # Use the same map coordinates as the input, with the reference latitude fixed ex ante.
    mask = available(pair)
    lat, x, y = pair["lat_deg"], pair["east_m"], pair["north_m"]
    candidates = mask & (abs(lat[:, None] + 9.25) <= 0.075)
    eligible_x = np.broadcast_to(x, mask.shape)[candidates]
    if len(eligible_x) == 0:
        return None, None
    median = np.median(eligible_x)
    nearest_lat = int(np.argmin(np.where(mask.any(axis=1), abs(lat + 9.25), np.inf)))
    center_x = x[np.argmin(np.where(mask[nearest_lat], abs(x - median), np.inf))]
    center_y = float(np.interp(-9.25, lat, y))
    distance2 = (x[None, :] - center_x)**2 + (y[:, None] - center_y)**2
    envelope = np.exp(-4 * np.log(2) * distance2 / 40000**2)
    return envelope, {"x_m": float(center_x), "lat_deg": -9.25, "fwhm_km": 40}


def synthetic_checks():
    lat = np.arange(-10.60, -7.90, 0.01)
    earlier = np.ones((2, len(lat), 18))
    later = np.ones_like(earlier)
    in_search = interval(lat, "search")
    later[:, in_search] *= 1.25
    pair = {"earlier": earlier, "later": later, "lat_deg": lat}
    assert profile(pair)["passes_gate"]
    assert not profile({**pair, "later": earlier})["passes_gate"]
    later[:, interval(lat, "north")] = np.nan
    assert not profile(pair)["passes_gate"]
    return ["injected_two_band_wide_brightening", "constant_scene_null", "missing_shoulder"]


def run(source):
    result = {"synthetic_checks": synthetic_checks(), "pairs": {}}
    with np.load(source, allow_pickle=False) as bundle:
        for name in ("primary", "control"):
            pair = load_pair(bundle, name, 720)
            stats = {
                "main_720m_linear": profile(pair),
                "sensitivity_720m_constant": profile(pair, constant=True),
                "sensitivity_360m_linear": profile(load_pair(bundle, name, 360)),
                "sensitivity_1440m_linear": profile(load_pair(bundle, name, 1440)),
            }
            passed = all(analysis["passes_gate"] for analysis in stats.values())
            run_data = {"maps": stats, "primary_and_sensitivity_pass": passed}
            run_data["registration_sensitivity"] = {
                direction: profile(shift(pair, dimension, delta))
                for direction, dimension, delta in (
                    ("east", 2, 1), ("west", 2, -1), ("north", 1, 1), ("south", 1, -1)
                )
            }
            if name == "primary":
                envelope, center = fake_patch(pair)
                run_data["injected_patch"] = center
                run_data["injected_results"] = {}
                if envelope is not None:
                    for amplitude in (0.05, 0.10, 0.20, 0.40):
                        injected = {**pair, "later": pair["later"] * (1 + amplitude * envelope)}
                        run_data["injected_results"][str(amplitude)] = profile(injected)
            result["pairs"][name] = run_data
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=ROOT / "matched-maps.npz")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = json.dumps(run(args.inputs), indent=2, allow_nan=False) + "\n"
    args.output.write_text(result, encoding="utf-8", newline="\n")
