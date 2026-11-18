"""Compare fixed Alice line contrasts under local and pooled dark subtraction."""

from argparse import ArgumentParser
from collections import defaultdict
import json
from pathlib import Path

import numpy as np


ROWS = {
    "stem": tuple(range(6, 16)),
    "stem_low": tuple(range(6, 11)),
    "stem_high": tuple(range(11, 16)),
    "edge_low": tuple(range(4)),
    "edge_high": tuple(range(28, 32)),
}
LINES = {
    "830": (830, 9, 1),
    "1135": (1134.6, 9, 1),
    "1135_wide": (1134.6, 15, 1),
    "1135_curved": (1134.6, 9, 2),
    "1549": (1549.5, 9, 1),
}


def bracket(observations):
    """Interpolate in time between separate closed exposures in each field."""
    grouped = defaultdict(list)
    for index, record in enumerate(observations):
        grouped[record["field"]].append(index)
    pairs, unbracketed, darks = {}, [], []
    for indices in grouped.values():
        field_darks = [index for index in indices if observations[index]["door"] == 1]
        darks.extend(field_darks)
        for index in indices:
            row = observations[index]
            if row["door"] != 2:
                continue
            compatible = [it for it in field_darks if same_configuration(observations[it], row)]
            before = [it for it in compatible if observations[it]["stop"] <= row["start"]]
            after = [it for it in compatible if observations[it]["start"] >= row["stop"]]
            if not before or not after:
                unbracketed.append(index)
                continue
            left, right = before[-1], after[0]
            t0, t1 = midpoint(observations[left]), midpoint(observations[right])
            right_weight = (midpoint(row) - t0) / (t1 - t0)
            pairs[index] = ((left, 1 - right_weight), (right, right_weight))
    return pairs, unbracketed, darks


def midpoint(row):
    return (row["start"] + row["stop"]) / 2


def same_configuration(a, b):
    return all(a[name] == b[name] for name in ("duration", "hv_safety", "hv_mask"))


def row_weights(wave, rows, center, width, degree):
    """Coefficients of a unit-height Gaussian relative to a local baseline."""
    contrast = np.zeros(wave.shape)
    for row in rows:
        x = (wave[row].astype(float) - center) / width
        support = (np.abs(x) <= 0.8) | ((np.abs(x) >= 1.6) & (np.abs(x) <= 4))
        selected = x[support]
        if min((selected < -1.6).sum(), (selected > 1.6).sum(), (np.abs(selected) <= 0.8).sum()) < 3:
            raise ValueError(f"Missing line or continuum coverage in detector row {row} at {center} Å")
        continuum = np.polynomial.polynomial.polyvander(selected, degree)
        gaussian = np.exp(-4 * np.log(2) * selected**2)
        design = np.column_stack((continuum, gaussian))
        contrast[row, support] = np.linalg.pinv(design)[-1] / len(rows)
    return contrast


def line_profiles(wave, rows, line):
    if line == "lya":
        profile = np.zeros(wave.shape)
        profile[:, rows, :] = ((wave[:, rows, :] >= 1195) & (wave[:, rows, :] <= 1235)) / len(rows)
        return profile
    return np.array([row_weights(frame, rows, *LINES[line]) for frame in wave])


def coefficients(profiles, observations, pairs, darks, sky):
    # Every paired sky exposure has the same duration by construction.
    n = len(sky)
    local = np.zeros(profiles.shape)
    pooled = np.zeros(profiles.shape)
    for index in sky:
        weights = profiles[index] / n
        duration = observations[index]["duration"]
        local[index] = weights / duration
        pooled[index] = local[index]
        for dark, fraction in pairs[index]:
            local[dark] -= fraction * weights / observations[dark]["duration"]
        compatible = [dark for dark in darks if same_configuration(observations[dark], observations[index])]
        exposure = sum(observations[dark]["duration"] for dark in compatible)
        # All exposures in the pool share this same acquisition configuration.
        pooled[compatible] -= weights / exposure
    return local, pooled


def measure(counts, coeff):
    value = float(np.einsum("ijk,ijk->", counts, coeff))
    variance = float(np.einsum("ijk,ijk,ijk->", counts, coeff, coeff))
    error = variance**0.5
    return {"rate": value, "count_error": error, "formal_ratio": value / error if error else None}


def audit(observations, counts, wave):
    pairs, unbracketed, darks = bracket(observations)
    subsets = {
        "all": sorted(pairs),
        "north": [i for i in pairs if observations[i]["hemisphere"] == "north"],
        "south": [i for i in pairs if observations[i]["hemisphere"] == "south"],
    }
    outputs = {}
    for aperture, rows in ROWS.items():
        outputs[aperture] = {}
        for line in (*LINES, "lya"):
            profiles = line_profiles(wave, rows, line)
            checks = {}
            for name, sky in subsets.items():
                local, pooled = coefficients(profiles, observations, pairs, darks, sky)
                checks[name] = {
                    "local": measure(counts, local),
                    "pooled": measure(counts, pooled),
                    "local_minus_pooled": measure(counts, local - pooled),
                }
            outputs[aperture][line] = checks
    profiles = line_profiles(wave, ROWS["stem"], "1135")
    field_results = {}
    for field in sorted({row["field"] for row in observations}):
        sky = [i for i in pairs if observations[i]["field"] == field]
        local, pooled = coefficients(profiles, observations, pairs, darks, sky)
        field_results[field] = {"sky": len(sky), "local": measure(counts, local),
                                "pooled": measure(counts, pooled)}
    field_rates = np.array([item["local"]["rate"] for item in field_results.values()])
    return {
        "sky_products": len(pairs) + len(unbracketed), "bracketed_sky_products": len(pairs),
        "dark_products": len(darks), "hemisphere_counts": {key: len(value) for key, value in subsets.items()},
        "unbracketed_mets": [observations[index]["met"] for index in unbracketed],
        "pairs": [
            {"met": observations[i]["met"],
             "darks": [dict(met=observations[j]["met"], weight=weight) for j, weight in pair]}
            for i, pair in pairs.items()
        ],
        "field_mean": float(field_rates.mean()),
        "field_mean_scatter_se": float(field_rates.std(ddof=1) / np.sqrt(len(field_rates))),
        "by_field": field_results,
        "contrasts": outputs,
    }


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    observations = json.loads((args.input_dir / "observations.json").read_bytes())["records"]
    with np.load(args.input_dir / "windows.npz", allow_pickle=False) as data:
        counts, wave = data["counts"], data["wave"]
    if counts.shape != wave.shape or counts.shape[:2] != (len(observations), 32):
        raise ValueError("Missing Alice detector data")
    results = audit(observations, counts, wave)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(results, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
