"""Audit the frozen all-slope mask and expand coverage with the slope held fixed."""

import argparse
import json
from pathlib import Path

import numpy as np

from analyze import MINUTES, SLOPES, change, delays, eligible, sample, summarize

ROOT = Path(__file__).resolve().parent


def run(data_path, summary_path):
    with np.load(data_path, allow_pickle=False) as inputs:
        solar, wind = inputs["solar"], inputs["wind"]
        velocity = inputs["flow"]
        distance = inputs["wind_position"] - inputs["solar_position"]
    selected = json.loads(summary_path.read_text(encoding="utf-8"))["selected_slope"]
    radial = int(np.flatnonzero(SLOPES == 0)[0]) + 1
    fitted = int(np.flatnonzero(SLOPES == selected)[0]) + 1
    models = [0, radial, fitted]

    normal = velocity[:, 0, None] + velocity[:, 1, None] * SLOPES
    finite_inputs = np.isfinite(distance).all(axis=1) & np.isfinite(velocity).all(axis=1)
    qualified = finite_inputs & (np.abs(normal) >= 100).all(axis=1)
    offset = delays(distance, velocity, SLOPES)
    geometry = np.isfinite(offset).all(axis=1)
    offset = np.column_stack((np.zeros(len(solar)), offset))
    matched = sample(wind, offset)
    solar_changes, wind_changes = change(solar), change(matched)
    all_models = eligible(solar_changes, wind_changes)
    three_models = eligible(solar_changes, wind_changes[:, models])
    solar_valid = np.isfinite(solar_changes).all(axis=1)
    both = np.zeros(len(solar), dtype=bool)
    both[5:] = geometry[5:] & geometry[:-5]
    both[np.arange(len(solar)) % MINUTES < 5] = False

    rows = []
    for day in np.split(np.arange(len(solar)), 3):
        rows.append({
            "qualified_flow_and_geometry_minutes": int(qualified[day].sum()),
            "qualified_minutes_failed_any_cap": int((qualified & ~geometry)[day].sum()),
            "geometry_eligible_changes": int(both[day].sum()),
            "excluded_missing_solar_changes": int((both & ~solar_valid)[day].sum()),
            "excluded_missing_wind_or_boundary": int(
                (both & solar_valid & ~all_models)[day].sum()
            ),
            "frozen_all_models": summarize(
                solar, matched, solar_changes, wind_changes, day, models, all_models
            ),
            "fixed_three_models": summarize(
                solar, matched, solar_changes, wind_changes, day, models, three_models
            ),
        })
    return {"selected_slope_from_frozen_summary": selected, "days": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=ROOT / "minute-inputs.npz")
    parser.add_argument("--summary", type=Path, default=ROOT / "summary.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    text = json.dumps(run(args.inputs, args.summary), allow_nan=False, indent=2) + "\n"
    args.output.write_text(text, encoding="utf-8", newline="\n")
