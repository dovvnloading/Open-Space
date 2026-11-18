"""Score Mariner 5 MAG boundary and cruise windows from the archived frames."""

import argparse
from bisect import bisect_left, bisect_right
from collections import Counter
import csv
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from hashlib import sha256
from math import ceil, isfinite, log, sqrt
from pathlib import Path
from statistics import median
import json
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
EPOCH = datetime(1967, 10, 19, 17, 34, 55, 300000, tzinfo=timezone.utc)
METRICS = ("magnitude", "vector_variability", "magnitude_variability")
FLOORS = (0.0, 0.2, 0.2)
REFERENCES = {"inbound": -156 * 60, "outbound": 18 * 60}


@dataclass(frozen=True)
class Frame:
    second: float
    identifier: int
    values: tuple


def relative_utc(value):
    time = datetime.fromisoformat(value.removesuffix("Z") + "+00:00")
    return (time - EPOCH).total_seconds()


def utc(second):
    return (EPOCH + timedelta(seconds=second)).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def encounter():
    path = ROOT / "mag_venus_encounter_data.csv"
    label = ET.parse(ROOT / "mag_venus_encounter_data.xml")
    declared = int(label.find(".//{*}Table_Delimited/{*}records").text)
    size = int(label.find(".//{*}file_size").text)
    receipt = next(r for r in json.loads((ROOT / "sources.json").read_bytes()) if r["file"] == path.name)
    body = path.read_bytes()
    if sha256(body).hexdigest() != receipt["raw_sha256"] or len(body) != size:
        raise ValueError("Encounter CSV does not match the receipt and label")
    with path.open(newline="", encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != declared:
        raise ValueError("Unexpected encounter record count")
    frames = []
    utc_errors = []
    excluded = 0
    for row in rows:
        second = float(row["HOUR"]) * 3600
        values = tuple(float(row[name]) for name in ("B", "SxASE", "SyASE", "SzASE", "SB"))
        if not isfinite(second) or any(not isfinite(n) or n < 0 for n in values) or values[0] <= 0:
            excluded += 1
            continue
        mag, dx, dy, dz, mag_std = values
        frames.append(Frame(second, int(row["FRAME_COUNT"]), (mag, sqrt(dx * dx + dy * dy + dz * dz), mag_std)))
        utc_errors.append(abs(relative_utc(row["UTC"]) - second))
    frames.sort(key=lambda f: f.second)
    frame_counts = Counter(f.second for f in frames)
    frames = [f for f in frames if frame_counts[f.second] == 1]
    return frames, {
        "declared_rows": declared,
        "retained_rows": len(frames),
        "invalid_rows": excluded,
        "duplicate_time_rows": sum(count for count in frame_counts.values() if count > 1),
        "largest_absolute_added_utc_difference_s": max(utc_errors),
        "gaps_over_80_s": sum(right.second - left.second > 80 for left, right in zip(frames, frames[1:])),
        "support_minutes": [frames[0].second / 60, frames[-1].second / 60],
    }


def coast():
    low, high = [], []
    with (ROOT / "cruise-controls.csv").open(newline="", encoding="ascii") as stream:
        for row in csv.DictReader(stream):
            second = relative_utc(row["utc"])
            values = tuple(float(row[key]) for key in ("mean_magnitude_nt", "vector_std_nt", "magnitude_std_nt"))
            frames = low if row["rate"] == "0" else high
            frames.append(Frame(second, int(row["record_n"]), values))
    return low, high


def sides(frames, times, at, near, far):
    before = frames[bisect_left(times, at - far):bisect_left(times, at - near)]
    after = frames[bisect_right(times, at + near):bisect_right(times, at + far)]
    return before, after


def measure(before, after, index, multiplier):
    left = median(frame.values[index] * multiplier for frame in before)
    right = median(frame.values[index] * multiplier for frame in after)
    score = log((right + FLOORS[index]) / (left + FLOORS[index]))
    return {
        "median_before_nt": left, "median_after_nt": right,
        "signed_log_ratio": score, "abs_log_ratio": abs(score),
        "count_before": len(before), "count_after": len(after),
    }


def scan(frames, center, near=120, far=480, minimum=5, limit_gap=80, deviation_scale=1):
    frames = sorted(frames, key=lambda f: f.second)
    times = [frame.second for frame in frames]
    maxima = dict.fromkeys(METRICS)
    floor = bisect_left(times, center - 600)
    ceiling = bisect_right(times, center + 600)
    valid_splits = 0
    for i in range(max(0, floor - 1), min(ceiling, len(frames) - 1)):
        left, right = frames[i:i + 2]
        at = (left.second + right.second) / 2
        if not center - 600 <= at <= center + 600 or not 0 < right.second - left.second <= limit_gap:
            continue
        before, after = sides(frames, times, at, near, far)
        if min(len(before), len(after)) < minimum:
            continue
        valid_splits += 1
        for index, key in enumerate(METRICS):
            result = measure(before, after, index, deviation_scale if index else 1)
            if maxima[key] is None or result["abs_log_ratio"] > maxima[key]["abs_log_ratio"]:
                maxima[key] = {
                    **result, "split_s": at, "split_utc": utc(at),
                    "offset_minutes": (at - center) / 60, "bracket_gap_s": right.second - left.second,
                    "left_frame": left.identifier, "right_frame": right.identifier,
                }
    before, after = sides(frames, times, center, near, far)
    fixed = None
    if min(len(before), len(after)) >= minimum:
        fixed = {key: measure(before, after, index, deviation_scale if index else 1)
                 for index, key in enumerate(METRICS)}
    return {"splits": valid_splits, "fixed_reference": fixed, "maxima": maxima}


def cruise_anchors():
    dates = [datetime(1967, 10, 1, tzinfo=timezone.utc) + timedelta(days=i) for i in range(15)]
    dates += [datetime(1967, 10, 22, tzinfo=timezone.utc) + timedelta(days=i) for i in range(15)]
    return [(day + timedelta(hours=hour) - EPOCH).total_seconds() for day in dates for hour in (3, 9, 15, 21)]


def scoring_context(controls):
    result = {}
    for metric in METRICS:
        values = sorted(row["scan"]["maxima"][metric]["abs_log_ratio"] for row in controls
                        if row["scan"]["maxima"][metric] is not None)
        result[metric] = {
            "eligible_blocks": len(values),
            "eligible_days": len({row["anchor_utc"][:10] for row in controls if row["scan"]["maxima"][metric]}),
            "p95_nearest_rank": values[ceil(0.95 * len(values)) - 1] if values else None,
        }
    return result


def compare(result, frames, controls, context, target):
    times = [frame.second for frame in frames]
    comparisons = {}
    for index, key in enumerate(METRICS):
        peak = result["maxima"][key]
        if peak is None:
            comparisons[key] = None
            continue
        expected_sign = 1 if target == "inbound" else -1
        threshold = context[key]["p95_nearest_rank"]
        value = peak["abs_log_ratio"]
        expected = peak["signed_log_ratio"] * expected_sign > 0
        within = abs(peak["offset_minutes"]) <= 5
        unusual = context[key]["eligible_blocks"] >= 60 and value > threshold
        before, after = sides(frames, times, peak["split_s"], 120, 480)
        deletions = []
        for side, other, is_before in ((before, after, True), (after, before, False)):
            for i in range(len(side)):
                candidate = side[:i] + side[i + 1:]
                if len(candidate) >= 5:
                    measured = measure(candidate, other, index, 1) if is_before else measure(other, candidate, index, 1)
                    deletions.append(measured)
        comparisons[key] = {
            "within_five_minutes": within,
            "in_expected_direction": expected,
            "exceeds_cruise_p95": unusual,
            "matched_either_sign": within and unusual,
            "matched_expected_sign": within and unusual and expected,
            "cruise_peaks_at_least_as_high": sum(
                row["scan"]["maxima"][key] is not None and row["scan"]["maxima"][key]["abs_log_ratio"] >= value
                for row in controls
            ),
            "deletion_possible": len(before) + len(after),
            "deletion_evaluable": len(deletions),
            "deletion_strength_and_original_sign": sum(
                threshold is not None and row["abs_log_ratio"] > threshold
                and row["signed_log_ratio"] * peak["signed_log_ratio"] > 0 for row in deletions
            ),
        }
    return comparisons


def run():
    frames, metadata = encounter()
    cruise, earth = coast()
    controls = [{"anchor_utc": utc(anchor), "scan": scan(cruise, anchor)} for anchor in cruise_anchors()]
    context = scoring_context(controls)
    events = {}
    for key, center in REFERENCES.items():
        primary = scan(frames, center)
        events[key] = {
            "reference_utc": utc(center),
            "primary": primary,
            "comparisons": compare(primary, frames, controls, context, key),
            "longer_windows": scan(frames, center, near=180, far=720, minimum=8),
            "sample_to_population_deviation_sensitivity": scan(frames, center, deviation_scale=sqrt(2 / 3)),
        }
    same_day = []
    for minute in range(-264, 49, 24):
        anchor = minute * 60
        if all(abs(anchor - center) > 28 * 60 for center in REFERENCES.values()):
            same_day.append({"anchor_utc": utc(anchor), "scan": scan(frames, anchor)})
    anchor = relative_utc("1967-06-14T12:23:00Z")
    return {
        "encounter": metadata,
        "events": events,
        "cruise_score_context": context,
        "cruise_anchors": controls,
        "encounter_day_controls": same_day,
        "earth_control": {"selected_frames": len(earth), "scan": scan(earth, anchor, minimum=20, limit_gap=25)},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
