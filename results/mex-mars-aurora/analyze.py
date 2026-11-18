"""Score Mars Express dual-band differences on geometry-defined limb tracks."""

import argparse
from bisect import bisect_left
from collections import defaultdict
import csv
from datetime import datetime
from io import StringIO
from pathlib import Path
from statistics import median, quantiles
import json
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
TARGET = (120, 210)
SKY = (350, 650)
TRIALS = ((SKY, "linear"), (SKY, "median"), ((300, 600), "linear"))


def second(value):
    return int(datetime.fromisoformat(value.rstrip("Z") + "+00:00").timestamp())


def load():
    with ZipFile(ROOT / "radio-seconds.zip") as archive:
        rows = list(csv.DictReader(StringIO(archive.read("radio-seconds.csv").decode("ascii"))))
    with ZipFile(ROOT / "geometry-2s.zip") as archive:
        frames = list(csv.DictReader(StringIO(archive.read("geometry-2s.csv").decode("ascii"))))
    geometry = defaultdict(dict)
    for pos in frames:
        geometry[pos["volume"]][second(pos["utc_reception"])] = tuple(float(pos[key]) for key in (
            "spherical_altitude_km", "planetocentric_lat_deg", "solar_zenith_deg", "ray_segment_fraction",
        ))
    data = defaultdict(list)
    for row in rows:
        start = second(row["utc"])
        coord = geometry[row["volume"]]
        # Each radio bin is [second, second+1), so geolocate it at the midpoint.
        even = start - start % 2
        if even not in coord or even + 2 not in coord:
            continue
        fraction = (start + 0.5 - even) / 2
        at, after = coord[even], coord[even + 2]
        altitude, latitude, sza, segment = (
            a + (b - a) * fraction for a, b in zip(at, after)
        )
        if not 0 <= segment <= 1:
            continue
        row.update({
            "second": start, "altitude": altitude, "latitude": latitude, "sza": sza,
            "direction": "ingress" if after[0] < at[0] else "egress",
            "d_hz": float(row["d_hz"]), "s_hz": float(row["s_hz"]), "x_hz": float(row["x_hz"]),
        })
        data[f"{row['volume']}-{row['chain']}"].append(row)
    return data


def groups(rows):
    if not rows:
        return []
    result = [[rows[0]]]
    for row in rows[1:]:
        previous = result[-1][-1]
        if row["second"] - previous["second"] > 2 or row["direction"] != previous["direction"]:
            result.append([])
        result[-1].append(row)
    return result


def inside(rows, interval):
    return [row for row in rows if interval[0] <= row["altitude"] < interval[1]]


def model(sample, column, method):
    origin = median(r["second"] for r in sample)
    values = [r[column] for r in sample]
    if method == "median":
        return origin, median(values), 0.0
    times = [r["second"] - origin for r in sample]
    center = sum(times) / len(times)
    average = sum(values) / len(values)
    slope = sum((t - center) * (y - average) for t, y in zip(times, values)) / sum((t - center) ** 2 for t in times)
    return origin, average - slope * center, slope


def distribute(values):
    return {"median_signed_hz": median(values),
            "median_absolute_hz": median(abs(x) for x in values),
            "p90_absolute_hz": quantiles([abs(x) for x in values], n=10, method="inclusive")[-1]}


def fit_branch(rows, background, method):
    sky, limb = inside(rows, background), inside(rows, TARGET)
    support = {"background_seconds": len(sky), "limb_seconds": len(limb)}
    if min(len(sky), len(limb)) < 15:
        return support
    r = {}
    for column in ("s_hz", "x_hz", "d_hz"):
        origin, y0, slope = model(sky, column, method)
        target = [row[column] - y0 - slope * (row["second"] - origin) for row in limb]
        controls = [row[column] - y0 - slope * (row["second"] - origin) for row in sky]
        ref = distribute(controls)
        limb_stats = distribute(target)
        ratio = limb_stats["median_absolute_hz"] / ref["median_absolute_hz"] if ref["median_absolute_hz"] else None
        r[column] = {"limb": limb_stats, "background": ref, "contrast": ratio, "slope_hz_per_s": slope}
    origin, y0, slope = model(sky, "d_hz", method)
    bins = defaultdict(list)
    for row in limb:
        bins[int(row["altitude"] // 10) * 10].append(row["d_hz"] - y0 - slope * (row["second"] - origin))
    return {
        **support, "median_latitude_deg": median(row["latitude"] for row in limb),
        "median_sza_deg": median(row["sza"] for row in limb),
        "limb_start_stop_utc": [limb[0]["utc"], limb[-1]["utc"]],
        "columns": r,
        "limb_10km_bins": {str(k): {"n": len(v), **distribute(v)} for k, v in sorted(bins.items()) if len(v) >= 3},
    }


def label_geometry(rows):
    errors = [abs(float(row["s_geometry_km"]) - (row["altitude"] + 3396.19))
              for row in inside(rows, (100, 700)) if row["s_geometry_km"]]
    return {"count": len(errors), "median_abs_difference_km": median(errors) if errors else None}


def pass_results(rows):
    results = []
    for group in groups(rows):
        if not inside(group, TARGET) and not inside(group, SKY):
            continue
        assessments = [fit_branch(group, background, approach) for background, approach in TRIALS]
        results.append({
            "direction": group[0]["direction"],
            "valid_seconds": len(group),
            "min_max_altitude_km": [min(r["altitude"] for r in group), max(r["altitude"] for r in group)],
            "primary": assessments[0], "offset_only": assessments[1], "alternate_sky": assessments[2],
        })
    return {
        "valid_geometry_matched_seconds": len(rows),
        "native_distance_comparison_near_limb": label_geometry(rows),
        "branches": results,
    }


def march_comparison(results):
    pairs = {}
    for chain in ("tcl1", "tcl2"):
        entries = []
        for volume in ("4432", "4434"):
            passes = results.get(f"{volume}-{chain}", {}).get("branches", [])
            entries.append(next((p["primary"] for p in passes
                                 if p["direction"] == "ingress" and "columns" in p["primary"]), None))
        before, after = entries
        if before is None or after is None:
            pairs[chain] = None
            continue
        sza = abs(after["median_sza_deg"] - before["median_sza_deg"])
        latitude = abs(after["median_latitude_deg"] - before["median_latitude_deg"])
        first, second_pass = (entry["limb_10km_bins"] for entry in entries)
        common = sorted(first.keys() & second_pass.keys(), key=int)
        difference = [second_pass[k]["median_absolute_hz"] - first[k]["median_absolute_hz"] for k in common]
        ratios = [entry["columns"]["d_hz"]["contrast"] for entry in entries]
        pairs[chain] = {
            "absolute_sza_difference_deg": sza, "absolute_latitude_difference_deg": latitude,
            "eligible_common_10km_bins": common,
            "geometry_gate_passes": len(common) >= 3 and sza <= 15 and latitude <= 20,
            "limb_to_background_contrast": ratios,
            "both_passes_twofold_contrast": all(x is not None and x >= 2 for x in ratios),
            "after_minus_before_abs_hz_by_bin": dict(zip(common, difference)),
            "median_after_minus_before_abs_hz": median(difference) if difference else None,
        }
    diffs = [pair["median_after_minus_before_abs_hz"] for pair in pairs.values()
             if pair is not None and pair["geometry_gate_passes"]]
    agreement = len(diffs) == 2 and diffs[0] * diffs[1] > 0
    return {"chains": pairs, "chains_agree_in_direction": agreement}


def run():
    results = {key: pass_results(rows) for key, rows in load().items()}
    return {"tracking_pairs": results, "march_17_18_new_norcia_comparison": march_comparison(results)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
