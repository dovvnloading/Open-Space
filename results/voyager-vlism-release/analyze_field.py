"""Score the frozen Voyager 1 plasma density and magnetic field comparison."""

import argparse
from collections import defaultdict
from datetime import date
import json
from math import log
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent
OFFSETS = ("-14", "-7", "7", "14")
BLOCK_ORIGIN = date(2023, 1, 1)
TRAIN_YEAR = "2023"
HOLDOUT_YEAR = "2025"


def day_values(rows, window, offset):
    days = defaultdict(list)
    for row in rows:
        field = row[window][offset]["field_nt"]
        if field is not None and row["low"] is not None and row["high"] is not None:
            days[row["scet"][:10]].append((row["density"], row["low"], row["high"], field))
    return {day: tuple(median(column) for column in zip(*values)) for day, values in days.items()}


def block_count(days):
    return len({(date.fromisoformat(day) - BLOCK_ORIGIN).days // 14 for day in days})


def score(days, predict):
    errors = []
    in_bounds = 0
    predictions = {}
    for day, (density, low, high, field) in sorted(days.items()):
        forecast = predict(field)
        errors.append(abs(log(forecast / density)))
        in_bounds += low <= forecast <= high
        predictions[day] = forecast
    return {
        "days": len(days),
        "blocks": block_count(days),
        "median_absolute_log_error": median(errors) if errors else None,
        "predictions_in_reported_bounds": in_bounds,
        "prediction_cm3_by_day": predictions,
    }


def subset_year(days, year):
    return {day: values for day, values in days.items() if day.startswith(year)}


def informative(result):
    return result["days"] >= 4 and result["blocks"] >= 3


def fit_window(rows, window):
    days = day_values(rows, window, "0")
    train = subset_year(days, TRAIN_YEAR)
    scale = median(density / field for density, _, _, field in train.values())
    baseline = median(density for density, _, _, _ in train.values())

    def predict(field):
        return scale * field

    yearly = {}
    for year in (TRAIN_YEAR, "2024", HOLDOUT_YEAR):
        sample = subset_year(days, year)
        yearly[year] = {
            "field_fit": score(sample, predict),
            "stationary_baseline": score(sample, lambda _: baseline),
            "simultaneous_by_day": {
                day: dict(zip(("density_cm3", "low_cm3", "high_cm3", "field_nt"), values))
                for day, values in sorted(sample.items())
            },
        }

    shifted = {}
    for offset in OFFSETS:
        shared_rows = [
            row for row in rows if row["scet"].startswith(HOLDOUT_YEAR)
            and row[window]["0"]["field_nt"] is not None
            and row[window][offset]["field_nt"] is not None
        ]
        shifted[offset] = {
            "simultaneous": score(day_values(shared_rows, window, "0"), predict),
            "shifted": score(day_values(shared_rows, window, offset), predict),
        }
    fit = yearly[HOLDOUT_YEAR]["field_fit"]
    reference = yearly[HOLDOUT_YEAR]["stationary_baseline"]
    beats_baseline = informative(fit) and fit["median_absolute_log_error"] <= (
        0.8 * reference["median_absolute_log_error"]
    )
    beats_all_offsets = all(
        informative(control["simultaneous"])
        and control["simultaneous"]["median_absolute_log_error"] < control["shifted"]["median_absolute_log_error"]
        for control in shifted.values()
    )
    return {
        "training_density_per_nt": scale,
        "training_stationary_density_cm3": baseline,
        "yearly": yearly,
        "held_out_shifted_controls_days": shifted,
        "simultaneous_holdout_informative": informative(fit),
        "beats_stationary_baseline_by_at_least_20_percent": beats_baseline,
        "beats_all_informative_shared_coverage_offsets": beats_all_offsets,
        "supports_frozen_model_rule": beats_baseline and beats_all_offsets,
    }


def embedded_control(rows):
    keys = ("-0.5", "0", "0.5")
    both = [
        row for row in rows if row["embedded"] is not None
        and all(row["embedded_w12"][key]["field_nt"] is not None for key in keys)
    ]
    return {
        "embedded_field_rows": sum(row["embedded"] is not None for row in rows),
        "shared_rows": len(both),
        "shared_utc_days": len({row["scet"][:10] for row in both}),
        "shared_median_absolute_difference_nt": {
            key: median(abs(row["embedded"] - row["embedded_w12"][key]["field_nt"]) for row in both)
            if both else None for key in keys
        },
    }


def run(data):
    rows = data["rows"]
    qtn = [row for row in rows if row["source"] == "qtn"]
    epo = [row for row in rows if row["source"] == "epo"]
    return {
        "observing_unit": "UTC day, with separate PWS rows retained in field-inputs.json",
        "block_anchor_utc": BLOCK_ORIGIN.isoformat(),
        "coverage": {
            year: {
                method: {
                    "rows": sum(row["scet"].startswith(year) for row in group),
                    "matched_rows_w12": sum(
                        row["scet"].startswith(year) and row["w12"]["0"]["field_nt"] is not None for row in group
                    ),
                } for method, group in (("qtn", qtn), ("epo", epo))
            } for year in (TRAIN_YEAR, "2024", HOLDOUT_YEAR)
        },
        "w12": fit_window(qtn, "w12"),
        "w6": fit_window(qtn, "w6"),
        "embedded_field_control": {
            "all": embedded_control(rows),
            "2025": embedded_control([row for row in rows if row["scet"].startswith(HOLDOUT_YEAR)]),
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "field-inputs.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    args.output.write_text(json.dumps(run(data), indent=2, allow_nan=False) + "\n",
                           encoding="utf-8", newline="\n")
