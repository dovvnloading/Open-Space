"""Evaluate frozen rest-wavelength predictions on a hash-pinned holdout roster."""

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import sys

from . import spectra
from .dataio import digest, load_observation
from .followup import measure
from .recurrence import C_KMS, _target_groups
from .windows import screen, window


def measure_spectrum(protocol, metadata, rows):
    redshift = metadata["redshift"]
    baseline = screen(rows, redshift)
    strict = screen([{**r, "mask": int(r["mask"]) | int(r["or_mask"])} for r in rows], redshift)
    tolerance = protocol["velocity_half_window_kms"]
    output = []
    for prediction in protocol["predictions"]:
        wave = prediction["vacuum_rest_angstrom"] * (1 + redshift)
        direction = prediction["direction"]
        primary = window(*baseline, wave, direction, tolerance)
        if primary["outcome"] == "hit":
            request = {"check": "continuum_sensitivity", "feature_id": primary["feature_id"]}
            fits = measure([request], baseline[2], rows)[0]["fits"]
            alternate = fits[1]
            primary["wide_continuum"] = alternate
            pixel = primary["center_pixel"]
            primary["core_or_mask"] = [int(row["or_mask"]) for row in rows[pixel - 1:pixel + 2]]
        controls = {}
        for offset in protocol["side_control_offsets_kms"]:
            shifted = wave * math.exp(offset / C_KMS)
            controls[str(offset)] = window(*baseline, shifted, direction, tolerance)
        output.append({"id": prediction["id"], "primary": primary, "shifted": controls,
                       "strict_mask": window(*strict, wave, direction, tolerance)})
    return output


def counts(by_observation):
    tally = Counter(row["outcome"] for _, row in by_observation)
    visits = sorted({label for label, row in by_observation if row["outcome"] == "hit"})
    return {"hits": tally["hit"], "covered": tally["hit"] + tally["miss"],
            "unsearchable": tally["unsearchable"], "visits_with_hits": visits}


def summarize(protocol, observations):
    summaries = []
    for i, prediction in enumerate(protocol["predictions"]):
        pairs = [(r["visit"], r["measurements"][i]) for r in observations]
        primary = counts([(visit, result["primary"]) for visit, result in pairs])
        strict = counts([(visit, result["strict_mask"]) for visit, result in pairs])
        offsets = {str(k): counts([(visit, r["shifted"][str(k)]) for visit, r in pairs])
                   for k in protocol["side_control_offsets_kms"]}
        fits = [r["primary"]["wide_continuum"] for _, r in pairs if r["primary"]["outcome"] == "hit"]
        contrasts = [fit["core_contrast"] for fit in fits if fit["status"] == "measured"]
        sign = 1 if prediction["direction"] == "emission" else -1
        summaries.append({**prediction, "primary": primary, "strict_mask": strict, "shifted": offsets,
                          "wide_continuum_evaluable_hits": len(contrasts),
                          "wide_continuum_same_sign_hits": sum(sign * c > 0 for c in contrasts),
                          "meets_recurrence_rule": (primary["hits"] >= protocol["minimum_targets"]
                                                    and len(primary["visits_with_hits"]) >= protocol["minimum_visits"])})
    return summaries


def run(protocol_path, roster_path, lock_path, input_directory):
    frozen = protocol_path.read_bytes()
    protocol = json.loads(frozen)
    if protocol["schema_version"] != 1 or digest(Path(spectra.__file__).read_bytes()) != protocol["detector_sha256"]:
        raise ValueError("Protocol schema or detector digest differs from the frozen experiment")
    lock = json.loads(lock_path.read_bytes())
    if (lock["protocol_sha256"] != digest(frozen)
            or lock["roster_sha256"] != digest(roster_path.read_bytes())):
        raise ValueError("Holdout protocol or roster differs from the sample lock")
    with roster_path.open(newline="") as stream:
        roster = list(csv.DictReader(stream))
    if not roster or len(roster) != lock["target_count"]:
        raise ValueError("Unexpected holdout roster size")
    targets = [{"ra_deg": r["ra"], "dec_deg": r["dec"]} for r in roster]
    if len(set(_target_groups(targets, 2.0))) != len(targets):
        raise ValueError("Repeated sky position in the holdout roster")
    observations = []
    ids = set()
    for record in roster:
        plate, mjd, fiber = (int(record[k]) for k in ("plate", "mjd", "fiber"))
        obs_id = f"sdss-{plate:04d}-{mjd:05d}-{fiber:04d}"
        if obs_id in ids or plate == protocol["discovery_visit"][0] or mjd == protocol["discovery_visit"][1]:
            raise ValueError(f"Repeated observation or discovery plate/night in holdout: {obs_id}")
        ids.add(obs_id)
        sidecar = input_directory / (obs_id + ".json")
        if digest(sidecar.read_bytes()) != record["manifest_sha256"]:
            raise ValueError(f"Observation metadata differs from roster for {obs_id}")
        meta, rows = load_observation(sidecar)
        if (meta["observation_id"] != obs_id or meta["source_sha256"] != record["source_sha256"]
                or meta["csv_sha256"] != record["csv_sha256"]
                or (meta["plate"], meta["mjd"], meta["fiber"]) != (plate, mjd, fiber)):
            raise ValueError(f"Observation identity or digest differs from roster for {obs_id}")
        observations.append({"observation_id": obs_id, "visit": f"{plate}-{mjd}", "redshift": meta["redshift"],
                             "csv_sha256": meta["csv_sha256"], "source_sha256": meta["source_sha256"],
                             "measurements": measure_spectrum(protocol, meta, rows)})
    return {"schema_version": 1, "claim_status": "screening_only",
            "protocol_sha256": digest(frozen), "roster_sha256": digest(roster_path.read_bytes()),
            "detector_sha256": protocol["detector_sha256"], "observation_count": len(observations),
            "summary": summarize(protocol, observations), "observations": observations}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, help="Directory of normalized observations")
    parser.add_argument("--protocol", type=Path, default=Path("protocols/holdout-v1.json"))
    parser.add_argument("--roster", type=Path, default=Path("data/holdout-targets.csv"))
    parser.add_argument("--lock", type=Path, default=Path("protocols/holdout-sample.json"))
    args = parser.parse_args(argv)
    try:
        result = run(args.protocol, args.roster, args.lock, args.inputs)
        print(json.dumps(result, indent=2, allow_nan=False))
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"holdout: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
