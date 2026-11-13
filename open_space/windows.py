"""Measure coverage and match signed features at predicted observed wavelengths."""

from bisect import bisect_left
import math

from . import spectra
from .recurrence import C_KMS


def screen(rows, redshift):
    pixels = spectra._pixels(rows)
    centers = {i for start, end in spectra._segments(pixels)
               for i in range(start + spectra.HALF_WINDOW, end - spectra.HALF_WINDOW)}
    analysis = spectra.analyze(rows, redshift)
    assert len(centers) == analysis["quality"]["searched_centers"]
    return [p[0] for p in pixels], centers, analysis["features"]


def window(grid, centers, features, predicted, direction, tolerance):
    upper = bisect_left(grid, predicted)
    nearest = min((i for i in (upper - 1, upper) if 0 <= i < len(grid)),
                  key=lambda i: (abs(predicted - grid[i]), i))
    result = {"predicted_observed_angstrom": predicted, "nearest_pixel": nearest}
    if nearest not in centers or not grid[0] <= predicted <= grid[-1]:
        return {**result, "outcome": "unsearchable"}
    matches = []
    for i, feature in enumerate(features):
        delta = C_KMS * math.log(feature["wavelength"] / predicted)
        if feature["direction"] == direction and abs(delta) <= tolerance:
            matches.append((abs(delta), -feature["score"], feature["center_pixel"], i, delta))
    if not matches:
        return {**result, "outcome": "miss"}
    _, _, center, index, delta = min(matches)
    hit = features[index]
    return {**result, "outcome": "hit", "feature_id": f"F{index + 1:04d}",
            "center_pixel": center, "observed_angstrom": hit["wavelength"], "offset_kms": delta,
            "score": hit["score"], "flux_contrast": hit["flux_contrast"],
            "catalog_classification": hit["classification"]}
