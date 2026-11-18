"""Replay the frozen SOLAR-1/Wind timing test from checked minute inputs."""

import argparse
import json
from pathlib import Path

import numpy as np


MINUTES = 1440
SLOPES = np.arange(-15, 16) / 10
ROOT = Path(__file__).resolve().parent


def delays(distance, velocity, slopes):
    normal = velocity[:, 0, None] + velocity[:, 1, None] * slopes
    travel = np.full_like(normal, np.nan)
    np.divide(distance[:, 0, None] + distance[:, 1, None] * slopes, normal,
              out=travel, where=np.abs(normal) >= 100)
    travel[np.abs(travel) > 3600] = np.nan
    return travel


def sample(wind, seconds):
    shifted = np.arange(len(wind), dtype=float)[:, None] + seconds / 60
    coordinate = np.where(np.isfinite(shifted), shifted, -2)
    nearest = np.rint(coordinate)
    coordinate = np.where(np.abs(coordinate - nearest) < 1e-10, nearest, coordinate)
    left = np.floor(coordinate).astype(int)
    right = np.ceil(coordinate).astype(int)
    native_day = np.arange(len(wind))[:, None] // MINUTES
    supported = (left // MINUTES == native_day) & (right // MINUTES == native_day)
    left = np.clip(left, 0, len(wind) - 1)
    right = np.clip(right, 0, len(wind) - 1)
    fraction = (coordinate - left)[..., None]
    result = (1 - fraction) * wind[left] + fraction * wind[right]
    result[~supported] = np.nan
    return result


def change(vectors):
    result = np.full_like(vectors, np.nan)
    result[5:] = vectors[5:] - vectors[:-5]
    result[np.arange(len(vectors)) % MINUTES < 5] = np.nan
    return result


def eligible(reference, comparisons):
    axes = tuple(range(1, comparisons.ndim))
    return np.isfinite(reference).all(axis=1) & np.isfinite(comparisons).all(axis=axes)


def score(reference, other):
    index = reference - np.mean(reference, axis=0)
    comparison = other - np.mean(other, axis=0)
    scale = np.linalg.norm(index) * np.linalg.norm(comparison)
    return float(np.vdot(index, comparison) / scale) if scale else None


def summarize(solar, wind, solar_change, wind_change, indices, models, keep, minimum=2):
    indices = indices[keep[indices]]
    if len(indices) < minimum:
        return {"n": len(indices), "score": None, "raw_vector_rms_nT": None}
    reference, raw = solar_change[indices], solar[indices]
    return {
        "n": len(indices),
        "score": [score(reference, wind_change[indices, model]) for model in models],
        "raw_vector_rms_nT": [
            float(np.sqrt(np.mean(np.sum((raw - wind[indices, model]) ** 2, axis=1))))
            for model in models
        ],
    }


def angle(first, second):
    norm = np.linalg.norm(first) * np.linalg.norm(second)
    if norm == 0:
        return None
    return float(np.rad2deg(np.arccos(np.clip(first @ second / norm, -1, 1))))


def rotation(magnetic, index):
    before = magnetic[index - 7:index - 1]
    after = magnetic[index + 2:index + 8]
    before = before[np.isfinite(before).all(axis=1)]
    after = after[np.isfinite(after).all(axis=1)]
    if min(len(before), len(after)) < 4:
        return None
    pre, post = np.median(before, axis=0), np.median(after, axis=0)
    if min(np.linalg.norm(pre), np.linalg.norm(post)) < 2:
        return None
    return angle(pre, post), post - pre


def describe_rotation(solar_delta, wind, index):
    if not 7 <= index < MINUTES - 7:
        return None
    signal = rotation(wind, index)
    if signal is None:
        return None
    amplitude, delta = signal
    disagreement = angle(solar_delta, delta)
    return {
        "minute": int(index),
        "utc": f"{index // 60:02d}:{index % 60:02d}:30",
        "rotation_deg": amplitude,
        "change_nT": float(np.linalg.norm(delta)),
        "change_direction_cosine": float(np.cos(np.deg2rad(disagreement)))
        if disagreement is not None else None,
    }


def rotations(solar, wind, model_delays):
    candidates = []
    for index in range(7, MINUTES - 7):
        estimate = rotation(solar, index)
        if estimate and estimate[0] >= 35:
            candidates.append((estimate[0], index, estimate[1]))
    selected = []
    for amplitude, index, delta in sorted(candidates, key=lambda row: (-row[0], row[1])):
        if any(abs(index - row["solar"]["minute"]) < 60 for row in selected):
            continue
        row = {"solar": describe_rotation(delta, solar, index), "wind": []}
        for travel in model_delays[index]:
            if not np.isfinite(travel):
                row["wind"].append(None)
                continue
            predicted = index + int(np.rint(travel / 60))
            direct = describe_rotation(delta, wind, predicted)
            nearby = [describe_rotation(delta, wind, minute)
                      for minute in range(predicted - 10, predicted + 11)]
            scored = [item for item in nearby
                      if item and item["change_direction_cosine"] is not None]
            peak = max(scored, key=lambda item: item["change_direction_cosine"], default=None)
            row["wind"].append({
                "predicted_minute": predicted, "direct": direct, "best_within_10_min": peak
            })
        selected.append(row)
        if len(selected) == 5:
            break
    return selected


def synthetic_checks():
    displacement = np.tile([-1800., 4200., 0.], (MINUTES, 1))
    velocity = np.tile([-600., 0., 0.], (MINUTES, 1))
    travel = delays(displacement, velocity, np.array([0., 1.]))
    assert np.allclose(travel[10], [3., -4.])
    rng = np.random.default_rng(1861)
    wind = rng.normal(size=(MINUTES, 3))
    planted = sample(wind, travel[:, 1:2])[:, 0]
    assert np.allclose(planted[10], wind[9] / 15 + wind[10] * 14 / 15)
    wind[9] = np.nan
    assert np.isnan(sample(wind, travel[:, 1:2])[10, 0]).all()
    assert np.array_equal(sample(wind, np.zeros((MINUTES, 1)))[10, 0], wind[10])
    assert np.isnan(sample(wind, travel[:, 0:1])[-1, 0]).all()
    assert score(change(planted)[20:-10], change(planted)[20:-10]) > 0.999999
    return ["radial_sign", "fractional_negative_lag", "mask_gap", "exact_center",
            "no_cross_day", "identical_change_score"]


def run(path):
    synthetic = synthetic_checks()
    with np.load(path, allow_pickle=False) as data:
        solar, wind = data["solar"], data["wind"]
        separation = data["wind_position"] - data["solar_position"]
        flow = data["flow"]
        counts = {
            key: [
                int(np.isfinite(series).all(axis=1).sum())
                for series in np.split(data[key], 3)
            ]
            for key in ("solar", "wind", "flow", "wind_position", "solar_position")
        }
        plasma_counts = data["plasma_counts"].tolist()
        gse_check = data["wind_gse_check_km"].tolist()
    travel = np.column_stack((np.zeros(len(solar)), delays(separation, flow, SLOPES)))
    grid = sample(wind, travel)
    solar_change, grid_change = change(solar), change(grid)
    shared = eligible(solar_change, grid_change)
    radial = int(np.flatnonzero(SLOPES == 0)[0]) + 1
    indices = np.arange(len(solar))
    training = indices[:2 * MINUTES]
    test = indices[2 * MINUTES:]
    model_grid = np.arange(travel.shape[1])
    fitted_training = summarize(solar, grid, solar_change, grid_change, training,
                                model_grid, shared)["score"]
    fit = max(range(1, len(fitted_training)), key=lambda model: (
        fitted_training[model], -abs(SLOPES[model - 1]), -SLOPES[model - 1]
    ))
    models = [0, radial, fit]
    daily = [
        summarize(solar, grid, solar_change, grid_change, day, models, shared)
        for day in np.split(indices, 3)
    ]
    blocks = [
        summarize(solar, grid, solar_change, grid_change, block, models, shared, minimum=80)
        for block in np.split(test, 12)
    ]
    qualifying = [block for block in blocks if block["score"] is not None]
    wins = sum(block["score"][2] > max(block["score"][:2]) for block in qualifying)
    heldout_score = daily[-1]["score"]
    controls = []
    for offset in range(120, 1321, 30):
        control = sample(np.roll(wind[test], offset, axis=0), travel[test][:, models])
        control_changes = change(control)
        cover = eligible(solar_change[test], control_changes) & shared[test]
        summary = summarize(solar[test], control, solar_change[test], control_changes,
                            np.arange(MINUTES), [0, 1, 2], cover, minimum=600)
        controls.append({"offset_min": offset, **summary})
    best_control = max(value for control in controls if control["score"] is not None
                       for value in control["score"])
    fixed = np.tile([-450., 0., 0.], (len(solar), 1))
    fixed_travel = delays(separation, fixed, SLOPES[[fit - 1]])
    fixed_grid = sample(wind, fixed_travel)
    fixed_changes = change(fixed_grid)
    fixed_keep = shared & eligible(solar_change, fixed_changes)
    fixed_result = summarize(solar, fixed_grid, solar_change, fixed_changes,
                             test, [0], fixed_keep)
    fixed_matched = summarize(solar, grid, solar_change, grid_change,
                              test, models, fixed_keep)
    same_subset_orbit = separation[test][shared[test]]
    heldout_delay = travel[test][shared[test]][:, models]
    return {
        "models": ["simultaneous", "radial", "fitted"],
        "coverage_valid_minutes": counts,
        "plasma_qc_counts": plasma_counts,
        "wind_gse_frame_check_max_abs_km": gse_check,
        "synthetic_checks": synthetic,
        "slopes": SLOPES.tolist(),
        "training_score_by_slope": fitted_training[1:],
        "selected_slope": float(SLOPES[fit - 1]),
        "days": ["2026-11-05", "2026-11-06", "2026-11-07"],
        "daily": daily,
        "heldout_blocks_2h": blocks,
        "heldout_qualifying_blocks": len(qualifying),
        "heldout_fitted_wins": wins,
        "heldout_delay_median_sec": np.median(heldout_delay, axis=0).tolist(),
        "heldout_delay_min_sec": np.min(heldout_delay, axis=0).tolist(),
        "heldout_delay_max_sec": np.max(heldout_delay, axis=0).tolist(),
        "heldout_wind_minus_solar_gse_median_km": np.median(same_subset_orbit, axis=0).tolist(),
        "heldout_transverse_separation_median_km": float(
            np.median(np.linalg.norm(same_subset_orbit[:, 1:], axis=1))
        ),
        "heldout_wind_vx_median_km_s": float(np.median(flow[test][shared[test], 0])),
        "heldout_fixed_speed_fitted": fixed_result,
        "heldout_fixed_speed_same_subset": fixed_matched,
        "heldout_circular_shift": controls,
        "heldout_peak_control_score": best_control,
        "heldout_improvement_gate": bool(
            heldout_score[2] - max(heldout_score[:2]) >= 0.05
            and wins > len(qualifying) / 2 and heldout_score[2] > best_control
        ),
        "heldout_rotations": rotations(solar[test], wind[test], travel[test][:, models]),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=ROOT / "minute-inputs.npz")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    output = json.dumps(run(arguments.inputs), allow_nan=False, indent=2) + "\n"
    arguments.output.write_text(output, encoding="utf-8", newline="\n")
