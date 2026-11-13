"""First-pass screening of narrow features on approximately regular optical grids.

Wavelengths are observed vacuum Angstroms. The dimensionless score orders
follow-up work. Its sampling distribution on survey data is uncalibrated.
"""

import math
import random
from statistics import median


# Source-frame rest wavelengths in vacuum Angstroms from the complete Table 2
# of the SDSS DR17 MaNGA pipeline, a limited set of optical/near-IR lines.
# https://www.sdss4.org/dr17/manga/manga-analysis-pipeline/
SOURCE_LINES = (
    ("[O II]", 3727.092), ("[O II]", 3729.875), ("H12", 3751.2174),
    ("H11", 3771.7012), ("H theta", 3798.9757), ("H eta", 3836.4720),
    ("[Ne III]", 3869.86), ("He I", 3889.749), ("H zeta", 3890.1506),
    ("[Ne III]", 3968.59), ("H epsilon", 3971.1951), ("H delta", 4102.8922),
    ("H gamma", 4341.6837), ("He II", 4687.015), ("H beta", 4862.683),
    ("[O III]", 4960.295), ("[O III]", 5008.24), ("[N I]", 5199.349),
    ("[N I]", 5201.705), ("He I", 5877.252), ("[O I]", 6302.046),
    ("[O I]", 6365.536), ("[N II]", 6549.86), ("H alpha", 6564.608),
    ("[N II]", 6585.27), ("[S II]", 6718.295), ("[S II]", 6732.674),
    ("He I", 7067.657), ("[Ar III]", 7137.76), ("[Ar III]", 7753.24),
    ("P eta", 9017.384), ("[S III]", 9071.10), ("P zeta", 9231.5460),
    ("[S III]", 9533.2), ("P epsilon", 9548.588),
)
# Three common airglow lines described in the SDSS DR17 optical caveats.
# https://www.sdss4.org/dr17/spectro/caveats/
SKY_LINES = (("[O I] airglow", 5578.9), ("[O I] airglow", 6302.046),
             ("[O I] airglow", 6365.536))
SCORE_THRESHOLD = 6.0
HALF_WINDOW = 12
SIDEBAND_START = 5
CATALOG_VERSION = "sdss-manga-dr17-table2-airglow3-v1"

# center, signed peak amplitude in units of the fixture's noise, pixel width
_INJECTIONS = ((4862.683, 15.0, 1.2), (5312.4, -12.0, 1.4),
               (5578.9, 17.0, 1.1), (5488.0, 60.0, 0.9))


def _fixture(seed, inject):
    rng = random.Random(seed)
    rows = []
    for wavelength in range(4800, 5641):
        continuum = 20.0 + 0.002 * (wavelength - 5200)
        flux = continuum + rng.gauss(0.0, 1.0)
        if inject:
            for center, amplitude, width in _INJECTIONS:
                flux += amplitude * math.exp(-0.5 * ((wavelength - center) / width) ** 2)
        rows.append({"wavelength": float(wavelength), "flux": flux, "ivar": 1.0,
                     "mask": int(abs(wavelength - 5488) <= 4)})
    return rows


def synthetic_rows(seed=0):
    """Return a reproducible injected fixture at z=0. No pixels are observations."""
    return _fixture(seed, inject=True)


def _pixels(rows):
    """Keep unusable samples in position so that masks break the search window."""
    pixels = []
    for index, row in enumerate(rows):
        try:
            wave = float(row["wavelength"])
            flux = float(row["flux"])
            weight = float(row["ivar"])
            # Accept integers and decimal integer strings from CSV. Fractional
            # masks are invalid, since truncation could silently clear a flag.
            raw_mask = row["mask"]
            mask = int(raw_mask)
            if not isinstance(raw_mask, str) and mask != raw_mask:
                raise ValueError("fractional mask")
        except (KeyError, TypeError, OverflowError, ValueError) as exc:
            raise ValueError(f"invalid spectrum fields at row {index + 1}: {exc}") from exc
        if not math.isfinite(wave) or wave <= 0 or (pixels and wave <= pixels[-1][0]):
            raise ValueError(f"wavelength must be positive, finite and increasing at row {index + 1}")
        usable = mask == 0 and math.isfinite(flux) and math.isfinite(weight) and weight > 0
        pixels.append((wave, flux, weight, usable))
    return pixels


def _segments(pixels):
    """Partition at every unusable pixel or spacing above 1.8 times the local median."""
    steps = [b[0] - a[0] for a, b in zip(pixels, pixels[1:])]
    start = 0
    for i, pixel in enumerate(pixels):
        if not pixel[3]:
            if i > start:
                yield start, i
            start = i + 1
        elif i > start:
            nearby = steps[max(0, i - 6):i - 1] + steps[i:min(len(steps), i + 5)]
            if nearby and steps[i - 1] > 1.8 * median(nearby):
                yield start, i
                start = i
    if start < len(pixels):
        yield start, len(pixels)


def _measure(pixels, i):
    """Compare the central triplet with a local continuum through two median anchors."""
    left = pixels[i - HALF_WINDOW:i - SIDEBAND_START + 1]
    right = pixels[i + SIDEBAND_START:i + HALF_WINDOW + 1]
    anchors = []
    for side in (left, right):
        anchors.append((median(p[0] for p in side), median(p[1] for p in side),
                        median(1.0 / p[2] for p in side) * math.pi / (2 * len(side))))
    (xl, yl, var_l), (xr, yr, var_r) = anchors

    def residual(p):
        return p[1] - (yl + (yr - yl) * ((p[0] - xl) / (xr - xl)))

    # Empirical scatter can inflate the model errors when the sidebands have
    # structure or the pipeline understates noise. MAD is robust to isolated peaks.
    normalized = [residual(p) * math.sqrt(p[2]) for p in left + right]
    offset = median(normalized)
    scale = max(1.0, 1.4826 * median(abs(n - offset) for n in normalized))
    core = pixels[i - 1:i + 2]
    total_weight = sum(p[2] for p in core)
    weighted_wave = sum(p[0] * (p[2] / total_weight) for p in core)
    residuals = [residual(p) for p in core]
    contrast = sum(r * (p[2] / total_weight) for r, p in zip(residuals, core))
    t = (weighted_wave - xl) / (xr - xl)
    variance = 1.0 / total_weight + (1 - t) ** 2 * var_l + t ** 2 * var_r
    score = abs(contrast) / (math.sqrt(variance) * scale)
    if not math.isfinite(score) or score < SCORE_THRESHOLD:
        return None
    sign = 1 if contrast > 0 else -1
    supported = [r * sign * math.sqrt(p[2]) / scale >= 2.0 for r, p in zip(residuals, core)]
    if not supported[1] or not (supported[0] or supported[2]):
        return None
    amplitudes = [max(0.0, r * sign) for r in residuals]
    area = sum(amplitudes)
    centroid = sum(p[0] * (a / area) for p, a in zip(core, amplitudes))
    return {"wavelength": centroid, "direction": "emission" if sign > 0 else "absorption",
            "score": score, "flux_contrast": contrast,
            "pixel_spacing_angstrom": (core[2][0] - core[0][0]) / 2.0,
            "center_pixel": i}


def _identify(feature, redshift):
    wave = feature["wavelength"]
    # A fixed pilot tolerance of two local pixels or 100 km/s, whichever is larger.
    # Includes typical barycentric motion in the observer-frame sky check.
    tolerance = max(2.0 * feature["pixel_spacing_angstrom"], wave * 100.0 / 299792.458)
    matches = []
    for frame, lines, factor in (("observer", SKY_LINES, 1.0),
                                 ("source", SOURCE_LINES, 1.0 + redshift)):
        for name, reference in lines:
            delta = wave - reference * factor
            if abs(delta) <= tolerance:
                matches.append({"name": name, "frame": frame, "reference_angstrom": reference,
                                "observed_offset_angstrom": delta})
    # Preserve source alternatives even when the sky veto takes precedence.
    if any(m["frame"] == "observer" for m in matches):
        classification = "sky_match"
    elif matches:
        classification = "catalog_match"
    else:
        classification = "unidentified"
    return {**feature, "rest_wavelength": wave / (1.0 + redshift), "matches": matches,
            "match_tolerance_angstrom": tolerance, "classification": classification}


def analyze(rows, redshift=0.0):
    """Screen an ascending spectrum. Nonzero masks and invalid flux/weights split coverage.

    Provide a measured source redshift (including a justified zero). Each segment
    needs at least 25 consecutive usable samples. A 'no_searchable_centers'
    result records absent coverage separately from a search yielding zero features.
    """
    redshift = float(redshift)
    if not math.isfinite(redshift) or redshift <= -1.0:
        raise ValueError("redshift must be finite and greater than -1")
    pixels = _pixels(rows)
    features = []
    searched = 0
    for start, end in _segments(pixels):
        cluster = []
        for i in range(start + HALF_WINDOW, end - HALF_WINDOW):
            searched += 1
            hit = _measure(pixels, i)
            # Collapse adjacent overlapping windows of the same sign into one feature.
            if cluster and (i - cluster[-1]["center_pixel"] > 2
                            or (hit and hit["direction"] != cluster[-1]["direction"])):
                features.append(_identify(max(cluster, key=lambda f: f["score"]), redshift))
                cluster = []
            if hit:
                cluster.append(hit)
        if cluster:
            features.append(_identify(max(cluster, key=lambda f: f["score"]), redshift))
    counts = {key: 0 for key in ("catalog_match", "sky_match", "unidentified")}
    for feature in features:
        counts[feature["classification"]] += 1
    return {
        "schema_version": 1, "method": "linear-sideband-triplet-v1", "catalog": CATALOG_VERSION,
        "redshift": redshift, "score_threshold": SCORE_THRESHOLD,
        "quality": {"status": "searched" if searched else "no_searchable_centers",
                    "input_pixels": len(pixels), "usable_pixels": sum(p[3] for p in pixels),
                    "searched_centers": searched},
        "features": features, "counts": counts,
        "interpretation": "Unidentified means outside the bundled line subset at the supplied "
                          "redshift. Scores rank narrow residuals and have no survey-wide "
                          "false-discovery calibration.",
    }


def benchmark(seeds=12):
    """Measure fixed-injection recovery and false triggers on separate simulated null spectra."""
    if isinstance(seeds, bool) or not isinstance(seeds, int) or seeds < 1:
        raise ValueError("seeds must be a positive integer")
    expected = (("catalog_match", 4862.683, "emission"), ("unidentified", 5312.4, "absorption"),
                ("sky_match", 5578.9, "emission"))
    recovered = {label: 0 for label, _, _ in expected}
    null_features = null_centers = masked_leaks = spurious = 0
    errors = []
    for seed in range(seeds):
        hits = analyze(synthetic_rows(seed))["features"]
        assigned = set()
        for label, center, direction in expected:
            candidates = [(abs(hit["wavelength"] - center), i) for i, hit in enumerate(hits)
                          if hit["classification"] == label and hit["direction"] == direction
                          and abs(hit["wavelength"] - center) <= 2.0]
            if candidates:
                error, index = min(candidates)
                recovered[label] += 1
                assigned.add(index)
                errors.append(error)
        masked_leaks += sum(abs(hit["wavelength"] - 5488.0) <= 5.0 for hit in hits)
        spurious += len(hits) - len(assigned)
        null = analyze(_fixture(1_000_000 + seed, inject=False))
        null_features += len(null["features"])
        null_centers += null["quality"]["searched_centers"]
    return {"kind": "synthetic_injection_recovery", "seeds": seeds,
            "injections_per_class": seeds, "recovered": recovered,
            "median_centroid_error_angstrom": median(errors) if errors else None,
            "unassigned_features_in_injection_runs": spurious, "masked_artifact_features": masked_leaks,
            "null_features": null_features, "null_searched_centers": null_centers}
