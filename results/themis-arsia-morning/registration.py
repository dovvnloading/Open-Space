"""Compare terrain texture outside the Arsia search interval at fixed translations."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, shift

from analyze import interval, load_pair

ROOT = Path(__file__).resolve().parent
SCALES = (360, 720)


def terrain(image, scale, shoulder_rows):
    valid = np.isfinite(image).all(axis=0) & shoulder_rows[:, None]
    sigma = 1440 / scale
    weight = gaussian_filter(valid.astype(float), sigma, mode="constant")
    result = []
    for band in image:
        leveled = np.zeros_like(band)
        np.log(band, where=valid, out=leveled)
        local = gaussian_filter(leveled, sigma, mode="constant")
        smooth = np.divide(local, weight, out=np.zeros_like(local), where=weight > 0)
        result.append(np.where(valid & (weight >= .8), leveled - smooth, np.nan))
    return np.stack(result)


def pearson_by_band(first, second):
    center_a, center_b = first - first.mean(axis=1)[:, None], second - second.mean(axis=1)[:, None]
    scale = np.linalg.norm(center_a, axis=1) * np.linalg.norm(center_b, axis=1)
    if np.any(scale == 0):
        raise ValueError("A terrain comparison is constant")
    return (np.sum(center_a * center_b, axis=1) / scale).tolist()


def score_on_shoulders(pair, scale):
    grid = range(-720, 721, scale)
    outputs = {}
    for shoulder in ("south", "north"):
        rows = interval(pair["lat_deg"], shoulder)
        first = terrain(pair["earlier"], scale, rows)
        second = terrain(pair["later"], scale, rows)
        shifts = {(x, y): shift(second, (0, y // scale, x // scale), order=0,
                                 mode="constant", cval=np.nan, prefilter=False)
                  for x in grid for y in grid}
        cells = np.isfinite(first).all(axis=0)
        for moved in shifts.values():
            cells &= np.isfinite(moved).all(axis=0)
        count = int(cells.sum())
        records = []
        outcome = {"identical_cells": count, "rows": records, "accepted_shift_m": None}
        if count >= 250:
            for (east, north), moved in shifts.items():
                by_band = pearson_by_band(first[:, cells], moved[:, cells])
                records.append({
                    "east_north_m": [east, north], "correlations": by_band,
                    "mean_correlation": float(np.mean(by_band)),
                })
            ranked = sorted(records, key=lambda row: (
                -row["mean_correlation"],
                sum(x*x for x in row["east_north_m"]),
                row["east_north_m"],
            ))
            winner, runner_up = ranked[:2]
            base = next(row for row in records if row["east_north_m"] == [0, 0])
            gain = winner["mean_correlation"] - base["mean_correlation"]
            separation = winner["mean_correlation"] - runner_up["mean_correlation"]
            outcome.update({
                "highest_score": winner,
                "zero_shift_score": base["mean_correlation"],
                "gain_over_zero": gain,
                "margin_to_runner_up": separation,
            })
            if winner["mean_correlation"] >= .3 and gain >= .06 and separation >= .015:
                outcome["accepted_shift_m"] = winner["east_north_m"]
        outputs[shoulder] = outcome
    return outputs


def test_no_external_pixels():
    image = np.ones((2, 20, 20))
    rows = (np.arange(20) >= 4) & (np.arange(20) < 15)
    changed = image.copy()
    changed[:, ~rows] *= 2
    np.testing.assert_array_equal(terrain(image, 720, rows), terrain(changed, 720, rows))


def run(source):
    test_no_external_pixels()
    with np.load(source, allow_pickle=False) as bundle:
        scopes = {str(scale): score_on_shoulders(load_pair(bundle, "primary", scale), scale)
                  for scale in SCALES}
    selected = [row["accepted_shift_m"] for data in scopes.values() for row in data.values()]
    translation = selected[0] if selected[0] and all(item == selected[0] for item in selected) else None
    return {"terrain_mask_check": "Outside-shoulder changes leave the terrain array identical",
            "scales": scopes, "jointly_supported_translation_m": translation}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=ROOT / "matched-maps.npz")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.dumps(run(args.inputs), indent=2, allow_nan=False) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
