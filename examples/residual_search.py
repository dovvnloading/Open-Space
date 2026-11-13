"""Blind narrow-line search of SDSS flux minus the fitted spectral model."""

import argparse
from bisect import bisect_left, bisect_right
from collections import Counter
import csv
import io
import json
import math
from pathlib import Path
from statistics import median
import sys

from open_space import recurrence, spectra
from open_space.archive import archive_source
from open_space.dataio import digest, load_observation


ROOT = Path(__file__).resolve().parents[1]
TOLERANCE_KMS = 150.0


def residual_pixels(rows, coadd_model, loglam):
    """Check coordinate alignment before subtracting each fitted flux density."""
    residuals = []
    for row, model, log_wave in zip(rows, coadd_model, loglam, strict=True):
        if not math.isclose(float(row["wavelength"]), 10 ** float(log_wave), rel_tol=0., abs_tol=1e-7):
            raise ValueError("COADD model grid differs from normalized wavelengths")
        residuals.append({**row, "flux": float(row["flux"]) - float(model)})
    return residuals


def searchable_waves(rows):
    p = spectra._pixels(rows)
    return [p[i][0] for start, stop in spectra._segments(p)
            for i in range(start + spectra.HALF_WINDOW, stop - spectra.HALF_WINDOW)]


def covered(waves, wavelength):
    """Does the ±150 km/s predicted window contain a searchable pixel center?"""
    factor = math.exp(TOLERANCE_KMS / recurrence.C_KMS)
    return bisect_right(waves, wavelength * factor) > bisect_left(waves, wavelength / factor)


def load_inputs(inputs, row):
    from astropy.io import fits

    obs = tuple(int(row[k]) for k in ("plate", "mjd", "fiber"))
    key, _ = archive_source(*obs)
    manifest = inputs / f"sdss-{key}.json"
    if digest(manifest.read_bytes()) != row["manifest_sha256"]:
        raise ValueError(f"Manifest digest mismatch for {key}")
    meta, pixels = load_observation(manifest)
    raw = (inputs / "raw" / f"{key}.fits").read_bytes()
    if (meta["observation_id"] != f"sdss-{key}" or digest(raw) != row["source_sha256"]
            or any(meta[k] != row[k] for k in ("source_sha256", "csv_sha256"))):
        raise ValueError(f"Input differs from pinned data for {key}")
    with fits.open(io.BytesIO(raw)) as fits_file:
        coadd = fits_file["COADD"].data
        residual = residual_pixels(pixels, coadd["model"], coadd["loglam"])
        sky = [float(v) for v in coadd["sky"]]
    return meta, pixels, residual, sky


def feature_evidence(feature, rows, sky, strict_scan):
    i = feature["center_pixel"]
    side = [sky[j] for j in (*range(i - 12, i - 4), *range(i + 5, i + 13)) if math.isfinite(sky[j])]
    local_sky = sky[i] if math.isfinite(sky[i]) else None
    side_median = median(side) if side else None
    core_or = [int(r["or_mask"]) for r in rows[i - 1:i + 2]]
    strict_coverage, strict_features = strict_scan
    matched = [f for f in strict_features if f["direction"] == feature["direction"] and
               abs(math.log(f["wavelength"] / feature["wavelength"]) * recurrence.C_KMS) <= TOLERANCE_KMS]
    return {"center_pixel": i, "flux_contrast": feature["flux_contrast"], "core_or_masks": core_or,
            "sky_at_center": local_sky, "sky_sideband_median": side_median,
            "sky_sideband_samples": len(side),
            "strict_or_window_covered": covered(strict_coverage, feature["wavelength"]),
            "strict_or_same_sign_match": bool(matched)}


def group_record(group, rank, coverages, frame):
    rest = frame == "rest"
    median_wave = group["median_rest_angstrom"] if rest else group["median_observed_angstrom"]
    available = []
    for obs_id, (z, waves) in coverages.items():
        prediction = median_wave * (1 + z) if rest else median_wave
        if covered(waves, prediction):
            available.append(obs_id)
    observed = [m["observation_id"] for m in group["members"]]
    # Centroids can sit slightly away from searchable pixel centers. Preserve the
    # window coverage and member counts separately instead of silently redefining either.
    return {"rank": rank, "direction": group["direction"], "targets": group["target_count"],
            "visits": group["visits"], "median_wavelength_angstrom": median_wave,
            "median_score": group["median_score"], "observed_span_angstrom": group["observed_span_angstrom"],
            "group_span_kms": group["frame_span_kms"], "searchable_targets_at_median": len(available),
            "member_targets_searchable_at_median": len(set(available).intersection(observed)),
            "passes_triage": (group["target_count"] >= 3 and
                              len({v["plate"] for v in group["visits"]}) >= 2 and
                              len({v["mjd"] for v in group["visits"]}) >= 2),
            "members": [{k: m[k] for k in ("observation_id", "feature_id", "direction", "score",
                                            "wavelength", "rest_wavelength")} for m in group["members"]]}


def run(inputs, roster):
    raw_roster = roster.read_bytes()
    with io.StringIO(raw_roster.decode("utf-8"), newline="") as stream:
        records = list(csv.DictReader(stream))
    if not records:
        raise ValueError("The input roster is empty")
    locations = [{"ra_deg": float(r["ra"]), "dec_deg": float(r["dec"])} for r in records]
    if len(set(recurrence._target_groups(locations, 2.0))) != len(records):
        raise ValueError("Repeated sky targets")
    packets, observations = [], []
    coverages = {}
    modules = [Path(__file__), Path(spectra.__file__), Path(recurrence.__file__)]
    modules_digest = {p.name: digest(p.read_bytes()) for p in modules}
    combined_sha = digest(json.dumps(modules_digest, sort_keys=True).encode())
    counts = {"raw": Counter(), "residual": Counter()}
    for row in records:
        meta, pixels, residuals, sky = load_inputs(inputs, row)
        target, z = meta["observation_id"], meta["redshift"]
        obs = {k: meta[k] for k in ("observation_id", "source_sha256", "csv_sha256", "redshift", "plate", "mjd")}
        for label, source in (("raw", pixels), ("residual", residuals)):
            scan = spectra.analyze(source, z)
            counts[label].update(scan["counts"])
            obs[label] = {"counts": scan["counts"], "quality": scan["quality"]}
            if label == "residual":
                scan["source_context"] = {"kind": "observed", "observation_id": target}
                packets.append({"schema_version": 1, "claim_status": "screening_only", "state": "completed",
                                "software_sha256": combined_sha, "observation": meta, "stages": {"detection": scan}})
                waves = searchable_waves(source)
                if len(waves) != scan["quality"]["searched_centers"]:
                    raise ValueError("Coverage disagrees with residual detector")
                coverages[target] = (z, waves)
                strict_rows = [{**r, "mask": int(r["mask"]) | int(r["or_mask"])} for r in source]
                strict_scan = (searchable_waves(strict_rows), spectra.analyze(strict_rows, z)["features"])
                obs["unassigned_residuals"] = [
                    {"feature_id": f"F{i:04d}", **{k: f[k] for k in ("wavelength", "rest_wavelength",
                                                                    "direction", "score")},
                     **feature_evidence(f, source, sky, strict_scan)}
                    for i, f in enumerate(scan["features"], 1) if f["classification"] == "unidentified"
                ]
        observations.append(obs)
    grouped = recurrence.summarize(packets, TOLERANCE_KMS)
    result = {"experiment": "blind_full_spectrum_sdss_model_residual_v1", "claim_status": "screening_only",
              "roster_sha256": digest(raw_roster), "implementation_sha256": modules_digest,
              "tolerance_kms": TOLERANCE_KMS, "catalog": spectra.CATALOG_VERSION,
              "score_threshold": spectra.SCORE_THRESHOLD, "triage_minimum": {"targets": 3, "plates": 2, "nights": 2},
              "targets": grouped["distinct_targets"], "counts": counts, "observations": observations}
    for frame, source_key in (("rest", "rest_frame_groups"), ("observer", "observer_frame_control_groups")):
        result[source_key] = [group_record(g, i, coverages, frame) for i, g in enumerate(grouped[source_key], 1)]
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, help="Directory holding the verified spectra and raw/ FITS files")
    parser.add_argument("--roster", type=Path, default=ROOT / "data/holdout-targets.csv")
    args = parser.parse_args(argv)
    try:
        result = run(args.inputs, args.roster)
        print(json.dumps(result, indent=2, allow_nan=False))
    except (ValueError, KeyError, TypeError, OSError, ImportError) as exc:
        print(f"residual-search: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
