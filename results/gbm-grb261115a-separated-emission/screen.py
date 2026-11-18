"""Repeat the frozen five-detector CTIME screen for Fermi GBM bn261115924."""

from argparse import ArgumentParser
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import platform

import astropy
from astropy.io import fits
import numpy as np
import scipy
from scipy.stats import binom


DETECTORS = ("n9", "n6", "n7", "n4", "n5")
REGIONS = {
    "pre": (-32, -3),
    "late": (8, 40),
    "pre_control": (-106, -77),
    "late_control": (77, 109),
}
SIDEBANDS = {"before": (-72, -36), "after": (44, 72)}
WIDTHS = (0.256, 0.512, 1.024, 2.048)
CADENCES = (0.064, 0.256)
TIME_TOLERANCE = 0.0001
ALPHA = 0.01
TRIGGER = 816473458.442344


@dataclass
class Detector:
    start: np.ndarray
    stop: np.ndarray
    counts: np.ndarray
    exposure: np.ndarray
    valid: np.ndarray
    band_kev: list


def read_detector(path, detector):
    with fits.open(path, memmap=False) as hdus:
        spec = hdus["SPECTRUM"]
        header = spec.header
        if (header["DETNAM"] != f"NAI_{int(detector[1:]):02}"
                or header["TIMESYS"] != "TT"
                or header["TIMEUNIT"] != "s"
                or header["TRIGTIME"] != TRIGGER):
            raise ValueError(f"Unexpected detector or time reference in {path.name}")
        rows = spec.data
        start = rows["TIME"] - TRIGGER
        stop = rows["ENDTIME"] - TRIGGER
        width = stop - start
        if np.any(width < 0) or np.any(start[1:] < stop[:-1] - TIME_TOLERANCE):
            raise ValueError(f"Invalid time grid in {path.name}")

        exposure = rows["EXPOSURE"].astype(float)
        mask = np.zeros(len(start), dtype=bool)
        for gti in hdus["GTI"].data:
            mask |= (start >= gti["START"] - TRIGGER - TIME_TOLERANCE) & (
                stop <= gti["STOP"] - TRIGGER + TIME_TOLERANCE
            )
        raw = rows["COUNTS"][:, 2:5]
        valid = (mask & (width > 0) & (rows["QUALITY"] == 0) & (exposure >= 0.90 * width)
                 & (exposure <= width + TIME_TOLERANCE) & np.all(raw >= 0, axis=1))
        edges = hdus["EBOUNDS"].data
        if edges["CHANNEL"].tolist() != list(range(8)):
            raise ValueError(f"Unexpected channel numbers in {path.name}")
        band = [[float(edges["E_MIN"][i]), float(edges["E_MAX"][i])] for i in range(2, 5)]
    return Detector(start, stop, raw.sum(axis=1), exposure, valid, band)


def full_bins(start, stop, limits):
    low, high = limits
    return (start >= low) & (stop <= high)


def summary(detector, limits):
    selected = full_bins(detector.start, detector.stop, limits)
    keep = selected & detector.valid
    exposure = float(detector.exposure[keep].sum())
    count = int(detector.counts[keep].sum())
    return {
        "native_bins": int(selected.sum()),
        "accepted_bins": int(keep.sum()),
        "covered_s": float((detector.stop[keep] - detector.start[keep]).sum()),
        "exposure_s": exposure,
        "counts": count,
        "rate": count / exposure if exposure else None,
    }


def exact_tail(count, exposure, background):
    bg_count = background["counts"]
    bg_exposure = background["exposure_s"]
    return float(binom.sf(count - 1, count + bg_count, exposure / (exposure + bg_exposure)))


def comparisons(count, exposure, backgrounds):
    p_before = exact_tail(count, exposure, backgrounds["before"])
    p_after = exact_tail(count, exposure, backgrounds["after"])
    return p_before, p_after, max(p_before, p_after)


def windows(start, stop, valid, limits):
    eligible = np.flatnonzero(full_bins(start, stop, limits) & valid)
    duration = stop - start
    for i in eligible:
        if i == 0:
            continue
        for nominal_width in WIDTHS:
            for cadence in CADENCES:
                n = round(nominal_width / cadence)
                j = i + n
                if j >= len(start) or stop[j - 1] > limits[1] or not valid[i:j].all():
                    continue
                if np.any(np.abs(duration[i - 1:j + 1] - cadence) > TIME_TOLERANCE):
                    continue
                if np.any(np.abs(start[i:j + 1] - stop[i - 1:j]) > TIME_TOLERANCE):
                    continue
                yield i, j, nominal_width


def score_windows(data, reference, common, backgrounds):
    output = []
    for region, limits in REGIONS.items():
        for i, j, width in windows(reference.start, reference.stop, common, limits):
            row = {"region": region, "start_s": float(reference.start[i]),
                   "stop_s": float(reference.stop[j - 1]), "nominal_width_s": width,
                   "native_bins": j - i}
            for name, detector in data.items():
                count = int(detector.counts[i:j].sum())
                exposure = float(detector.exposure[i:j].sum())
                before, after, score = comparisons(count, exposure, backgrounds[name])
                row.update({
                    f"{name}_counts": count, f"{name}_exposure_s": exposure,
                    f"{name}_p_before": before, f"{name}_p_after": after, f"{name}_p": score,
                })
            row["joint_p"] = max(row["n9_p"], row["n6_p"])
            output.append(row)

    trials = sum(row["region"] in ("pre", "late") for row in output)
    if not trials:
        raise ValueError("No complete primary search windows")
    for row in output:
        row["primary_scaled_p"] = min(1.0, row["joint_p"] * trials)
        row["n7_passes"] = row["n7_p"] <= ALPHA
        row["direction_passes"] = min(row["n4_p"], row["n5_p"]) > ALPHA
        row["passes"] = (row["primary_scaled_p"] <= ALPHA
                         and row["n7_passes"] and row["direction_passes"])
    return output, trials


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(Path(__file__).with_name("inputs.json").read_text(encoding="utf-8"))
    for name, expected in manifest["files"].items():
        if hashlib.sha256((args.inputs / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Digest mismatch in {name}")

    data = {
        name: read_detector(args.inputs / f"glg_ctime_{name}_bn261115924_v00.pha", name)
        for name in DETECTORS
    }
    reference = data["n9"]
    for name, detector in data.items():
        if not (np.array_equal(detector.start, reference.start)
                and np.array_equal(detector.stop, reference.stop)):
            raise ValueError(f"Unaligned time grid in {name}")
    common = np.logical_and.reduce([detector.valid for detector in data.values()])

    backgrounds = {
        name: {side: summary(detector, interval) for side, interval in SIDEBANDS.items()}
        for name, detector in data.items()
    }
    for sides in backgrounds.values():
        for name, (low, high) in SIDEBANDS.items():
            if sides[name]["exposure_s"] < 0.90 * (high - low):
                raise ValueError(f"Insufficient exposure in {name} sideband")

    output, trials = score_windows(data, reference, common, backgrounds)
    region_results = {}
    for region, limits in REGIONS.items():
        selected = [row for row in output if row["region"] == region]
        eligible = full_bins(reference.start, reference.stop, limits)
        region_results[region] = {
            "limits_s": limits,
            "total_native_bins": int(eligible.sum()),
            "joint_valid_native_bins": int((eligible & common).sum()),
            "by_width": {str(width): sum(row["nominal_width_s"] == width for row in selected)
                         for width in WIDTHS},
            "trials": len(selected),
            "primary_threshold_crossings": sum(row["primary_scaled_p"] <= ALPHA for row in selected),
            "all_gate_crossings": sum(row["passes"] for row in selected),
            "best": min(selected, key=lambda row: row["joint_p"]) if selected else None,
        }

    prompt = {name: summary(detector, (-1, 4)) for name, detector in data.items()}
    subintervals = {}
    for name, detector in data.items():
        prompt[name]["p_before"], prompt[name]["p_after"], prompt[name]["p"] = comparisons(
            prompt[name]["counts"], prompt[name]["exposure_s"], backgrounds[name]
        )
        subintervals[name] = [
            {"limits_s": [int(left), int(right)], **summary(detector, (left, right))}
            for band in SIDEBANDS.values()
            for left, right in zip(np.linspace(*band, 5)[:-1], np.linspace(*band, 5)[1:])
        ]

    result = {
        "trigger_met": TRIGGER,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                    "astropy": astropy.__version__, "scipy": scipy.__version__},
        "family_trials": trials,
        "family_alpha": ALPHA,
        "unadjusted_joint_threshold": ALPHA / trials,
        "excluded_zero_duration_rows": {
            name: int((detector.stop == detector.start).sum()) for name, detector in data.items()
        },
        "energy_bounds_kev": {name: detector.band_kev for name, detector in data.items()},
        "backgrounds": backgrounds,
        "background_subintervals": subintervals,
        "prompt_control": prompt,
        "regions": region_results,
    }
    args.output.mkdir(parents=True)
    with (args.output / "windows.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    (args.output / "screen.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n"
    )


if __name__ == "__main__":
    main()
