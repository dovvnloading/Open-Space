"""Score Antarctic ozone threshold area and fixed-location persistence."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

ROOT = Path(__file__).resolve().parent
YEARS = (2024, 2025, 2026)
THRESHOLDS = (170, 180, 190, 200, 220, 240)
SENSORS = ("omi", "omps")
LAT_EDGES = np.deg2rad(np.arange(-90, -39))
WEIGHTS = np.broadcast_to(6371 ** 2 * np.deg2rad(1) * np.diff(np.sin(LAT_EDGES))[:, None], (50, 360))
CAP_70 = np.indices(WEIGHTS.shape)[0] < 20
CAP_60 = np.indices(WEIGHTS.shape)[0] < 30
CAP_40 = np.ones(WEIGHTS.shape, dtype=bool)


def area(mask):
    return float(WEIGHTS[mask].sum()) / 1e6


def read():
    metadata = json.loads((ROOT / "polar-metadata.json").read_bytes())
    data = {}
    with ZipFile(ROOT / "polar-grids.zip") as bundle:
        for member, item in metadata["members"].items():
            body = bundle.read(member)
            if sha256(body).hexdigest() != item["polar_sha256"]:
                raise ValueError(f"Polar member digest mismatch: {member}")
            date, sensor = member.removesuffix(".u16").rsplit("-", 1)
            data.setdefault(date, {})[sensor] = np.frombuffer(body, dtype="<u2").reshape(50, 360)
    return {date: np.stack([pair[sensor] for sensor in SENSORS]).astype(np.int32) for date, pair in sorted(data.items())}


def weighted_median(values, weights):
    if not len(values):
        return None
    order = np.argsort(values)
    return float(values[order][np.searchsorted(np.cumsum(weights[order]), weights.sum() / 2)])


def measure_day(values, region):
    valid = (values > 0) & region
    common = np.all(valid, axis=0)
    offsets = values[0] - values[1]
    near = common & (abs(values[0] + values[1] - 440) <= 30)
    outs = {}
    for threshold in (180, 200, 220, 240):
        low = valid & (values < threshold)
        established = np.all(low, axis=0)
        possible = np.all(low | ~valid, axis=0) & region
        disagrees = common & (low[0] != low[1])
        outs[str(threshold)] = {
            "common_omi_low_mkm2": area(common & low[0]),
            "common_omps_low_mkm2": area(common & low[1]),
            "omi_only_low_mkm2": area(common & low[0] & ~low[1]),
            "omps_only_low_mkm2": area(common & ~low[0] & low[1]),
            "union_on_common_mkm2": area(common & (low[0] | low[1])),
            "consensus_lower_mkm2": area(established),
            "consensus_upper_mkm2": area(possible),
            "omi_observed_low_mkm2": area(low[0]),
            "omps_observed_low_mkm2": area(low[1]),
            "common_agreed_low_cells": int(established.sum()),
        }
        if threshold == 220:
            outs[str(threshold)].update({
                "near_threshold_disagree_mkm2": area(near & disagrees),
                "near_threshold_common_mkm2": area(near),
                "omi_common_deficit_du_million_km2": float((WEIGHTS[common] * np.maximum(0, 220 - values[0][common])).sum() / 1e6),
                "omps_common_deficit_du_million_km2": float((WEIGHTS[common] * np.maximum(0, 220 - values[1][common])).sum() / 1e6),
            })
    return {
        "cap_mkm2": area(region),
        "omi_valid_mkm2": area(valid[0]),
        "omps_valid_mkm2": area(valid[1]),
        "common_mkm2": area(common),
        "common_cells": int(common.sum()),
        "median_omi_minus_omps_du": weighted_median(offsets[common], WEIGHTS[common]),
        "median_near_220_omi_minus_omps_du": weighted_median(offsets[near], WEIGHTS[near]),
        "thresholds": outs,
    }


def components(mask):
    candidates = set(zip(*np.nonzero(mask), strict=True))
    result = []
    while candidates:
        seed = candidates.pop()
        cells = [seed]
        for row, col in cells:
            for point in ((row - 1, col), (row + 1, col), (row, (col - 1) % 360), (row, (col + 1) % 360)):
                if point in candidates:
                    candidates.remove(point)
                    cells.append(point)
        weights = WEIGHTS[tuple(zip(*cells, strict=True))]
        result.append((float(weights.sum()) / 1e6, cells))
    return sorted(result, key=lambda x: (-x[0], min(x[1])))


def split_check(values):
    valid = values > 0
    low = valid & (values < 220)
    strict = components(np.all(low, axis=0))
    possible = np.all(low | ~valid, axis=0)
    stats = {
        "consensus_largest_two_mkm2": [group[0] for group in strict[:2]],
        "omi_largest_two_mkm2": [group[0] for group in components(low[0])[:2]],
        "omps_largest_two_mkm2": [group[0] for group in components(low[1])[:2]],
        "second_at_least_half_mkm2": len(strict) >= 2 and strict[1][0] >= 0.5,
        "separation_survives_possible_missing_bridge": None,
    }
    if stats["second_at_least_half_mkm2"]:
        first, second = [group[1][0] for group in strict[:2]]
        shares = any(first in group and second in group for _, group in components(possible))
        stats["separation_survives_possible_missing_bridge"] = not shares
    return stats


def support(windows):
    return [np.all(window > 0, axis=1) for window in windows]


def thresholds(windows, footprint, relaxed=False):
    early, late = windows
    _, late_valid = support(windows)
    summary = {}
    for threshold in THRESHOLDS:
        initial = footprint & ((np.all((early > 0) & (early < threshold), axis=1)).sum(axis=0) >= 4)
        late_low = np.all((late > 0) & (late < threshold), axis=1).sum(axis=0)
        enough = late_valid.sum(axis=0) if relaxed else 5
        final = footprint & (late_low == enough)
        initial_area, same_area = area(initial), area(initial & final)
        summary[str(threshold)] = {
            "early_area_mkm2": initial_area,
            "late_area_mkm2": area(final),
            "retained_area_mkm2": same_area,
            "early_only_mkm2": area(initial & ~final),
            "late_only_mkm2": area(final & ~initial),
            "fraction_of_early_area_retained": same_area / initial_area if initial_area else None,
        }
    return summary


def coverage(mask):
    return {
        "retained_mkm2": area(mask),
        "fraction_of_60s_cap": area(mask) / area(CAP_60),
        "by_five_degree_band_mkm2": [
            float(WEIGHTS[i:i + 5][mask[i:i + 5]].sum()) / 1e6 for i in range(0, 30, 5)
        ],
        "by_sixty_degree_sector_mkm2": [
            float(WEIGHTS[:, i:i + 60][mask[:, i:i + 60]].sum()) / 1e6 for i in range(0, 360, 60)
        ],
    }


def run(data):
    windows = {year: [
        np.stack([data[f"{year}-11-{day:02d}"] for day in range(first, first + 5)])
        for first in (1, 11)
    ] for year in YEARS}
    fixed = {year: CAP_60 & np.all(np.concatenate(support(sample)), axis=0) for year, sample in windows.items()}
    across = np.all(list(fixed.values()), axis=0)
    study = {}
    for year, sample in windows.items():
        early_valid, late_valid = support(sample)
        relaxed = CAP_60 & (early_valid.sum(axis=0) >= 4) & (late_valid.sum(axis=0) >= 4)
        study[str(year)] = {
            "within_year_footprint": coverage(fixed[year]),
            "relaxed_footprint": coverage(relaxed),
            "fixed_across_years": thresholds(sample, across),
            "fixed_within_year": thresholds(sample, fixed[year]),
            "fixed_70s_within_year": thresholds(sample, fixed[year] & CAP_70),
            "relaxed_within_year": thresholds(sample, relaxed, relaxed=True),
        }
    per_day = {date: {"south_60": measure_day(maps, CAP_60), "south_40": measure_day(maps, CAP_40)}
               for date, maps in data.items()}
    splits = {date: split_check(data[date]) for date in data if date.startswith("2026-11-")}
    qualified = [date for date, row in splits.items() if row["separation_survives_possible_missing_bridge"] is True]
    return {
        "across_years_footprint": coverage(across),
        "seasons": study,
        "daily": per_day,
        "november_2026_split_screen": splits,
        "consecutive_robust_split_pairs": [[a, b] for a, b in zip(qualified, qualified[1:])
                                             if int(b[-2:]) - int(a[-2:]) == 1],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(read()), indent=2, allow_nan=False) + "\n",
                           encoding="utf-8", newline="\n")
