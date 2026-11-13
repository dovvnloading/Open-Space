"""Paired artificial-line recovery experiments on an immutable observed spectrum."""

import argparse
import json
import math
import random
import sys
from pathlib import Path
from statistics import median

from .dataio import load_observation
from .spectra import CATALOG_VERSION, HALF_WINDOW, SCORE_THRESHOLD, _pixels, _segments, analyze


AMPLITUDES = (3.0, 6.0, 10.0)
SIGNS = (("emission", 1), ("absorption", -1))
MATCH_PIXELS = 2.0
SIGMA_RANGE_PIXELS = (0.9, 1.7)


def _inject(rows, pixels, index, wave, width, signed_peak):
    """Copy only rows touched by a Gaussian truncated at six profile sigma."""
    trial = list(rows)
    for step in (-1, 1):
        j = index if step == -1 else index + 1
        while 0 <= j < len(pixels):
            position, flux, _, usable = pixels[j]
            delta = (position - wave) / width
            if abs(delta) > 6.0:
                break
            if usable:
                trial[j] = {**rows[j], "flux": flux + signed_peak * math.exp(-0.5 * delta ** 2)}
            j += step
    return trial


def _recover(truth_wave, tolerance, direction, baseline, measured):
    def near(wave, reference):
        return abs(wave - reference) <= tolerance

    if any(near(hit["wavelength"], truth_wave) for hit in baseline):
        return "baseline_confounded", None
    matched = [hit["wavelength"] for hit in measured
               if hit["direction"] == direction and near(hit["wavelength"], truth_wave)]
    if not matched:
        return "missed", None
    recovered_wave = min(matched, key=lambda wave: abs(wave - truth_wave))
    if any(near(hit["wavelength"], recovered_wave) for hit in baseline):
        return "baseline_confounded", None
    return "recovered", recovered_wave - truth_wave


def calibrate(manifest, rows, trials=24, seed=17):
    """Sample positions with replacement, pairing amplitudes and signs at each position.

    Centers derive only from the grid, mask and variance. The detector operates
    blind to truth on the entire perturbed spectrum. The unchanged baseline enters
    only the post-detection accounting. All results belong to calibration.
    """
    if isinstance(trials, bool) or not isinstance(trials, int) or trials < 1:
        raise ValueError("trials must be a positive integer")
    rows = list(rows)
    pixels = _pixels(rows)
    centers = [i for start, end in _segments(pixels)
               for i in range(start + HALF_WINDOW, end - HALF_WINDOW)]
    if not centers:
        raise ValueError("the spectrum has no searchable injection centers")
    z = float(manifest["redshift"])
    baseline = analyze(rows, z)
    rng = random.Random(seed)
    outcomes = []
    results = {}
    for direction, _ in SIGNS:
        for amplitude in AMPLITUDES:
            results[direction, amplitude] = {
                "direction": direction, "peak_noise_units": amplitude, "attempted": trials,
                "recovered": 0, "missed": 0, "baseline_confounded": 0, "offsets": [],
            }
    for _ in range(trials):
        i = rng.choice(centers)
        central_wave, _, ivar, _ = pixels[i]
        spacing = (pixels[i + 1][0] - pixels[i - 1][0]) / 2.0
        wave = central_wave + rng.uniform(-0.45, 0.45) * spacing
        sigma_pixels = rng.uniform(*SIGMA_RANGE_PIXELS)
        position = {"pixel": i, "injected_wavelength": wave, "sigma_pixels": sigma_pixels,
                    "pixel_spacing_angstrom": spacing, "outcomes": []}
        for direction, sign in SIGNS:
            for amplitude in AMPLITUDES:
                perturbed = _inject(rows, pixels, i, wave, sigma_pixels * spacing,
                                    sign * amplitude / math.sqrt(ivar))
                measured = analyze(perturbed, z)
                status, offset = _recover(wave, MATCH_PIXELS * spacing, direction,
                                          baseline["features"], measured["features"])
                cell = results[direction, amplitude]
                cell[status] += 1
                if offset is not None:
                    cell["offsets"].append(abs(offset))
                position["outcomes"].append({"direction": direction, "peak_noise_units": amplitude,
                                             "outcome": status, "offset_angstrom": offset})
        outcomes.append(position)
    for cell in results.values():
        assessable = trials - cell["baseline_confounded"]
        cell["assessable"] = assessable
        cell["recovered_fraction_assessable"] = cell["recovered"] / assessable if assessable else None
        offsets = cell.pop("offsets")
        cell["median_absolute_offset_angstrom"] = median(offsets) if offsets else None
    return {
        "schema_version": 1, "experiment": "real_background_injection_recovery_v1",
        "data_role": "calibration", "observation_id": manifest["observation_id"],
        "source_sha256": manifest["source_sha256"], "csv_sha256": manifest["csv_sha256"],
        "redshift": z, "detector": baseline["method"], "catalog": CATALOG_VERSION,
        "score_threshold": SCORE_THRESHOLD, "seed": seed, "trials_per_cell": trials,
        "injected_spectra": trials * len(SIGNS) * len(AMPLITUDES),
        "searchable_centers": len(centers), "input_pixels": len(pixels),
        "usable_pixels": baseline["quality"]["usable_pixels"],
        "baseline_features": len(baseline["features"]),
        "distinct_center_pixels_sampled": len({r["pixel"] for r in outcomes}),
        "centroid_tolerance_pixels": MATCH_PIXELS, "summary": list(results.values()),
        "positions": outcomes,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="path to a digest-pinned observation manifest")
    parser.add_argument("--trials", type=int, default=24, help="positions per sign/amplitude cell")
    parser.add_argument("--seed", type=int, default=17, help="seed for position and width draws")
    args = parser.parse_args(argv)
    try:
        manifest, rows = load_observation(args.manifest)
        result = calibrate(manifest, rows, trials=args.trials, seed=args.seed)
    except (KeyError, OSError, ValueError) as exc:
        parser.exit(2, f"calibration failed: {exc}\n")
    json.dump(result, sys.stdout, sort_keys=True, allow_nan=False, indent=2)
    print()


if __name__ == "__main__":
    main()
