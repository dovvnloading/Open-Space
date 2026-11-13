"""Replay the fixed-position Ca II K and [O III] test on the original campaign's repeat visit."""

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import random
import re
from statistics import median, quantiles
from urllib.parse import urlencode
from urllib.request import urlopen

from examples.calcium_doublet import K
from open_space import archive
from open_space.cohort import CATALOG, read_bounded

ROOT = Path(__file__).resolve().parents[1]
CATALOG_HASH = "7dda43a170ed86b12d380b824babdbb756748a3677c194efd7cfcb6a1a42b5d3"
FITS_HASHES = {
    (6497, 56329, 139): archive.PILOT_SHA256,
    (6012, 56101, 585): "7903cca82ee86f1d17e4e2432ac27b808e5e67be0908b4c86ae60b75affae44e",
}
LINES = {"Ca_II_K": K, "K_minus_25": K - 25, "K_plus_25": K + 25, "O_III_5008": 5008.239}


def separation(ra1, dec1, ra2, dec2):
    """Great-circle distance in arcseconds, stable at coincident positions."""
    a, b = math.radians(dec1), math.radians(dec2)
    chord = math.sin((a-b)/2)**2 + math.cos(a)*math.cos(b)*math.sin(math.radians(ra1-ra2)/2)**2
    return math.degrees(2 * math.asin(math.sqrt(min(1, chord)))) * 3600


def paired_targets(raw):
    """Return alternate coadds for the original 12, keeping their anchor identifiers."""
    records = list(csv.DictReader(line for line in raw.decode().splitlines() if not line.startswith("#")))
    anchors = sorted({int(row["anchor"]) for row in records})
    expected = [int(item["observation_id"].split("-")[-1])
                for item in json.loads((ROOT / "results/first-visit.json").read_text())["observations"]]
    if anchors != sorted(expected):
        raise ValueError("Catalog missing one or more original observations")
    pairs = []
    for row in records:
        obs = (int(row["plate"]), int(row["mjd"]), int(row["fiber"]))
        if obs[1] == 56329 or int(row["bestobjid"]) <= 0:
            continue
        offset = separation(*(float(row[k]) for k in ("anchor_ra", "anchor_dec", "ra", "dec")))
        if (offset <= 2 and row["run2d"] == archive.RUN2D and row["survey"] == "boss"
                and row["class_noqso"] == "GALAXY" and int(row["zwarning_noqso"]) == 0):
            pairs.append(((6497, 56329, int(row["anchor"])), obs, row["bestobjid"], offset))
    return anchors, pairs


def exposure_ids(header):
    """Reduce blue/red camera records to science integrations, retaining the original cards."""
    cards = sorted((key, str(header[key]).strip()) for key in header if re.fullmatch(r"EXPID\d+", key))
    if not cards or len(cards) != header["NEXP"]:
        raise ValueError("Missing or incomplete exposure headers")
    ids = set()
    for _, value in cards:
        match = re.fullmatch(r"[br][12]-(\d{8})-\d{8}-\d{8}", value)
        if match is None:
            raise ValueError(f"Unrecognized exposure record {value}")
        ids.add(match[1])
    return {"science_ids": sorted(ids), "camera_records": dict(cards)}


def fixed_window(rows, z, center, mask_column, noise_seed=None):
    """Return signed absorption depth and coverage, using the anchor z for both visits."""
    window = []
    for row in rows:
        rest = float(row["wavelength"]) / (1 + z)
        distance = rest - center
        if abs(distance) <= 15:
            valid = int(row[mask_column]) == 0 and float(row["ivar"]) > 0
            window.append((distance, float(row["flux"]), float(row["ivar"]), valid))
    sides = [[p for p in window if start <= p[0] <= stop and p[3]]
             for start, stop in ((-15., -8.), (8., 15.))]
    core = [p for p in window if abs(p[0]) <= 3. and p[3]]
    record = {"pixels_blue_core_red": [len(sides[0]), len(core), len(sides[1])]}
    if min(len(s) for s in sides) < 2 or len(core) < 3:
        return {**record, "status": "insufficient_pixels"}
    x0, x1 = (median(p[0] for p in side) for side in sides)
    fractions = [(p[0] - x0) / (x1 - x0) for p in core]
    total_weight = sum(p[2] for p in core)

    def estimate(blue, central, red):
        y0, y1 = median(blue), median(red)
        continuum = [y0 + (y1 - y0) * f for f in fractions]
        if min(continuum) <= 0:
            return None
        mean_continuum = sum(c * p[2] for c, p in zip(continuum, core)) / total_weight
        mean_flux = sum(f * p[2] for f, p in zip(central, core)) / total_weight
        return mean_flux, mean_continuum, 1 - mean_flux / mean_continuum

    windows = (sides[0], core, sides[1])
    estimate0 = estimate(*([p[1] for p in group] for group in windows))
    if estimate0 is None:
        return {**record, "status": "nonpositive_continuum"}
    flux, continuum, depth = estimate0
    record.update(status="measured", weighted_flux=flux, weighted_continuum=continuum,
                  flux_contrast=flux - continuum, fractional_depth=depth)
    if noise_seed is not None:
        rng = random.Random(noise_seed)
        depths = []
        for _ in range(4096):
            simulated = ([rng.gauss(p[1], 1 / math.sqrt(p[2])) for p in group] for group in windows)
            estimate_i = estimate(*simulated)
            if estimate_i is not None:
                depths.append(estimate_i[2])
        if len(depths) >= 2:
            percentiles = quantiles(depths, n=1000, method="inclusive")
            record["diagonal_noise_sensitivity"] = {
                "draws": 4096, "usable": len(depths), "seed": noise_seed,
                "p025": percentiles[24], "p975": percentiles[974],
            }
    return record


def spectrum(obs, cache):
    from astropy.io import fits

    key, url = archive.archive_source(*obs)
    path = cache / f"spec-{key}.fits"
    raw, sha = archive.read_source(url, path if path.exists() else None, FITS_HASHES[obs])
    archive.write_immutable(path, raw)
    normalized, details = archive.csv_and_metadata(raw, obs)
    with fits.open(io.BytesIO(raw)) as hdus:
        exposure = exposure_ids(hdus[0].header)
        thing = int(hdus["SPALL"].data[0]["THING_ID"])
    if details["classification"] != "GALAXY" or details["zwarning"] != 0 or details["plate_quality"] != "good":
        raise ValueError("Repeat test requires good-plate, zero-warning galaxy solutions")
    provenance = {"observation_id": f"sdss-{key}", "source_url": url, "source_sha256": sha,
                  "csv_sha256": hashlib.sha256(normalized).hexdigest(), "thing_id": thing,
                  "plate": obs[0], "mjd": obs[1], "fiber": obs[2], **details, **exposure}
    return provenance, list(csv.DictReader(io.StringIO(normalized.decode())))


def run(cache, refresh_catalog=False):
    cache.mkdir(parents=True, exist_ok=True)
    sql = (ROOT / "data/repeat-catalog.sql").read_text()
    if refresh_catalog:
        with urlopen(CATALOG + "?" + urlencode({"cmd": sql, "format":"csv"}), timeout=60) as response:
            raw = read_bounded(response)
    else:
        raw = (ROOT / "data/repeat-catalog.csv").read_bytes()
    catalog_hash = hashlib.sha256(raw).hexdigest()
    if catalog_hash != CATALOG_HASH:
        raise ValueError("Catalog response differs from the frozen metadata selection")
    anchors, pairs = paired_targets(raw)
    if len(pairs) != 1:
        raise ValueError("Expected one eligible catalog pair in frozen data")
    anchor, repeat, photometric_id, catalog_offset = pairs[0]
    anchor_meta, anchor_rows = spectrum(anchor, cache)
    repeat_meta, repeat_rows = spectrum(repeat, cache)
    z = anchor_meta["redshift"]
    separation_arcsec = separation(anchor_meta["ra_deg"],anchor_meta["dec_deg"],repeat_meta["ra_deg"],repeat_meta["dec_deg"])
    if (anchor_meta["thing_id"] <= 0 or anchor_meta["thing_id"] != repeat_meta["thing_id"]
            or separation_arcsec > 2):
        raise ValueError("FITS positions or target identifiers differ")
    overlap = sorted(set(anchor_meta["science_ids"]) & set(repeat_meta["science_ids"]))
    for meta, rows in ((anchor_meta, anchor_rows), (repeat_meta, repeat_rows)):
        meta["indices"] = {}
        for name, rest in LINES.items():
            meta["indices"][name] = {
                "rest_vacuum_angstrom": rest, "predicted_observed_angstrom": rest * (1+z),
                "default": fixed_window(rows, z, rest, "and_mask", 17 if name == "Ca_II_K" else None),
                "strict": fixed_window(rows, z, rest, "or_mask", 17 if name == "Ca_II_K" else None),
            }
        oxygen = meta["indices"]["O_III_5008"]
        for mask, policy in (("and_mask","default"), ("or_mask","strict")):
            candidates = [r for r in rows if int(r[mask]) == 0 and float(r["ivar"]) > 0
                          and abs(float(r["wavelength"]) - oxygen["predicted_observed_angstrom"]) <= 4]
            eligible = max(candidates, key=lambda r: float(r["flux"])) if candidates else None
            oxygen[policy]["eligible_peak_pixel_angstrom"] = float(eligible["wavelength"]) if eligible else None
    return {"schema_version":1, "kind":"same_target_repeat", "catalog_url":CATALOG,
            "catalog_sha256":catalog_hash, "catalog_anchors":anchors, "eligible_repeat_pairs":len(pairs),
            "photometric_bestobjid":photometric_id, "catalog_offset_arcsec":catalog_offset,
            "fit_position_offset_arcsec":separation_arcsec, "days_between_coadd_mjd":abs(anchor[1]-repeat[1]),
            "reference_redshift":z, "shared_science_exposures":overlap,
            "exposure_status":"overlap" if overlap else "disjoint",
            "spectra":[anchor_meta,repeat_meta]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("runs/repeat-fits"))
    parser.add_argument("--refresh-catalog", action="store_true", help="Verify live SkyServer response")
    args = parser.parse_args()
    print(json.dumps(run(args.cache, args.refresh_catalog), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
