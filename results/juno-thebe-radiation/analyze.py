"""Replay the frozen two-detector Juno experiment from compact checked inputs."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
CONTROLS = (-20, -15, -10, 10, 15, 20)
PIXELS = ("T1_shielded", "T3_shielded", "T2_bare", "T4_bare", "T5_bare")


def binned_jedi(data, center, strict):
    time, duration, dead, mode, quality, sector = data[:, :6].T
    counts = data[:, 6:]
    offsets = time - center
    good = (
        (np.abs(offsets) < 2700) & (mode == 0) & (quality == 0)
        & np.isfinite(duration) & (duration > 0)
        & np.isfinite(counts).all(axis=1) & (counts >= 0).all(axis=1)
    )
    if strict:
        good &= (dead >= 0) & (dead <= 0.1 * duration)
    selected = np.flatnonzero(good)
    bins = np.floor(offsets[selected] / 30).astype(int)
    centers = np.arange(-90, 90, dtype=int) * 30 + 15
    rates = np.full((len(centers), len(PIXELS)), np.nan)
    support = []
    for index, offset in enumerate(centers):
        rows = selected[bins == (offset - 15) // 30]
        seconds = float(duration[rows].sum())
        sectors = len(np.unique(sector[rows]))
        qualified = seconds >= 5 and sectors >= 4
        if qualified:
            rates[index] = np.sum(counts[rows], axis=0) / seconds
        support.append({
            "offset_s": int(offset), "rows": len(rows), "duration_s": seconds,
            "sectors": sectors, "qualified": qualified,
        })
    return centers, rates, support, int(good.sum())


def local_contrast(centers, rates, offset):
    times = centers - offset
    valid = np.isfinite(rates).all(axis=1)
    central = rates[valid & (np.abs(times) < 60)]
    earlier = rates[valid & (times >= -1200) & (times <= -480)]
    later = rates[valid & (times >= 480) & (times <= 1200)]
    counts = [len(central), len(earlier), len(later)]
    result = {"bin_counts_center_before_after": counts, "relative_deficit": None}
    if np.any(np.array(counts) < (3, 8, 8)):
        return result
    pre = np.median(earlier, axis=0)
    post = np.median(later, axis=0)
    center = np.median(central, axis=0)
    baseline = np.sqrt(pre * post)
    result["before_cps"] = pre.tolist()
    result["after_cps"] = post.tolist()
    result["central_cps"] = center.tolist()
    result["baseline_cps"] = baseline.tolist()
    if np.any(baseline <= 0):
        return result
    result["relative_deficit"] = (1 - center / baseline).tolist()
    return result


def jedi_test(data, et, strict):
    centers, rates, support, retained = binned_jedi(data, et, strict)
    contrasts = {str(minute): local_contrast(centers, rates, minute * 60)
                 for minute in (0, *CONTROLS)}
    event = contrasts["0"]["relative_deficit"]
    controls = [contrasts[str(minute)]["relative_deficit"] for minute in CONTROLS]
    supported = [trial[:2] for trial in controls if trial is not None]
    peak = np.max(supported, axis=0) if supported else None
    passed = (
        event is not None and len(supported) >= 4
        and np.all(np.asarray(event[:2]) >= 0.30)
        and np.all(np.asarray(event[:2]) > peak)
    )
    for index, record in enumerate(support):
        record["cps"] = rates[index].tolist() if record["qualified"] else None
    return {
        "rows_in_range": len(data), "retained_rows": retained, "bins": support,
        "contrasts_by_minute_center": contrasts, "supported_controls": len(supported),
        "peak_shielded_control_deficit": peak.tolist() if peak is not None else None,
        "passes_gate": bool(passed),
    }


def describe_sru(orbit):
    sru = orbit["sru"]
    images = sru["images"]
    jupiter, thebe = orbit["jupiter"], orbit["thebe"]
    times = np.array([image["et"] for image in images])
    nearest_index = int(np.argmin(abs(times - thebe["et"])))
    nearest_jupiter = int(np.argmin(abs(times - jupiter["et"])))
    triplet = images[max(0, nearest_index - 1):nearest_index + 2]
    identical = (len(triplet) == 3 and len({image["sensor"] for image in triplet}) == 1
                 and len({image["method"] for image in triplet}) == 1
                 and all(image["derived"] is not None and image["method"] is not None
                         for image in triplet))
    deficit = None
    if identical:
        before, center, after = [image["derived"] for image in triplet]
        deficit = [1 - center / before, 1 - center / after]
    return {
        "table_rows": len(images), "label_rows": sru["declared_rows"],
        "thebe_nearest_image": images[nearest_index],
        "thebe_offset_sec": float(times[nearest_index] - thebe["et"]),
        "jupiter_nearest_image": images[nearest_jupiter],
        "jupiter_offset_sec": float(times[nearest_jupiter] - jupiter["et"]),
        "adjacent_exposure_window": triplet, "relative_deficit_vs_each_neighbor": deficit,
        "passes_matching_mode_gate": identical,
        "passes_deficit_gate": deficit is not None and min(deficit) >= .30,
    }


def compare_camera(asc, sru):
    time, counts, bbo = asc.T
    good = (bbo == 0) & np.isfinite(counts) & (counts > 0)
    matches = []
    for image in sru["images"]:
        instant = image["utc_seconds_since_j2000_noon"]
        within = good & (np.abs(time - instant) <= 120)
        covered = np.min(time) <= instant <= np.max(time) and within.sum() >= 10
        matches.append({
            "utc": image["utc"], "sru_derived_rate": image["derived"],
            "sru_method": image["method"], "asc_rows": int(within.sum()),
            "asc_median_calibrated_counts": float(np.median(counts[within])) if covered else None,
        })
    joint = [(row["sru_derived_rate"], row["asc_median_calibrated_counts"]) for row in matches
             if row["sru_derived_rate"] is not None and row["asc_median_calibrated_counts"] is not None]
    rank = None
    if len(joint) >= 5 and (np.ptp(joint, axis=0) > 0).all():
        rank = float(spearmanr(*np.asarray(joint).T).statistic)
    return {"matches": matches, "paired_rows": len(joint), "spearman_rho": rank}


def test_local_contrast():
    time = np.arange(-90, 90) * 30 + 15
    rates = np.full((len(time), len(PIXELS)), 10.)
    rates[np.abs(time) < 60, :2] = 5.
    expected = local_contrast(time, rates, 0)
    assert expected["relative_deficit"] == [0.5, 0.5, 0., 0., 0.]
    rates[(np.abs(time) < 60) & (time > 0)] = np.nan
    assert local_contrast(time, rates, 0)["relative_deficit"] is None
    return ["injected_half_rate_in_both_witnesses", "two_bin_central_coverage_failure"]


def run(bundle, source):
    geometry = json.loads(source.read_text(encoding="utf-8"))
    result = {"pixel_order": PIXELS, "synthetic_checks": test_local_contrast(), "passes": {}}
    with np.load(bundle, allow_pickle=False) as data:
        for orbit, row in geometry.items():
            watch = {"geometry": {key: row[key] for key in ("jupiter", "thebe")},
                     "sru": describe_sru(row)}
            if orbit in ("81", "82"):
                watch["sru_asc"] = compare_camera(data[f"asc_{orbit}"], row["sru"])
            if orbit in ("82", "83"):
                center = row["thebe" if orbit == "83" else "jupiter"]["et"]
                series = data[f"jedi_{orbit}"]
                watch["jedi_primary"] = jedi_test(series, center, strict=False)
                watch["jedi_low_deadtime"] = jedi_test(series, center, strict=True)
            if orbit == "83":
                watch["all_detector_gates_pass"] = bool(
                    watch["jedi_primary"]["passes_gate"] and watch["jedi_low_deadtime"]["passes_gate"]
                    and watch["sru"]["passes_deficit_gate"]
                )
            result["passes"][orbit] = watch
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=ROOT / "compact-inputs.npz")
    parser.add_argument("--geometry", type=Path, default=ROOT / "geometry-sru.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = json.dumps(run(args.inputs, args.geometry), allow_nan=False, indent=2) + "\n"
    args.output.write_text(report, encoding="utf-8", newline="\n")
