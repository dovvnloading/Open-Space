"""Explore how much detected local contrast survives subtracting the SDSS fit."""

import argparse
import json
import math
from pathlib import Path
from statistics import median
import sys

from open_space.dataio import digest, load_observation
from open_space.followup import measure
from open_space.windows import screen, window


def core_contrast(rows, pixel):
    fit = measure([{"feature_id": "F0001", "check": "continuum_sensitivity"}],
                  [{"center_pixel": pixel, "wavelength": float(rows[pixel]["wavelength"])}], rows)[0]["fits"][0]
    return fit.get("core_contrast")


def compare(evidence, protocol, inputs):
    from astropy.io import fits

    results = {p["id"]: [] for p in protocol["predictions"]}
    for obs in evidence["observations"]:
        obs_id = obs["observation_id"]
        meta, rows = load_observation(inputs / f"{obs_id}.json")
        raw_path = inputs / "raw" / (obs_id.removeprefix("sdss-") + ".fits")
        if (meta["source_sha256"] != obs["source_sha256"] or digest(raw_path.read_bytes()) != obs["source_sha256"]
                or meta["csv_sha256"] != obs["csv_sha256"] or meta["redshift"] != obs["redshift"]):
            raise ValueError(f"Input differs from original measurement for {obs_id}")
        with fits.open(raw_path) as hdus:
            model = [float(x) for x in hdus["COADD"].data["model"]]
        model_rows = [{**r, "flux": v} for r, v in zip(rows, model, strict=True)]
        residual_rows = [{**r, "flux": float(r["flux"]) - v} for r, v in zip(rows, model, strict=True)]
        residual_screen = screen(residual_rows, meta["redshift"])
        for p, measured in zip(protocol["predictions"], obs["measurements"], strict=True):
            primary = measured["primary"]
            if p["id"] != measured["id"]:
                raise ValueError("Prediction order changed")
            if primary["outcome"] != "hit":
                continue
            pixel = primary["center_pixel"]
            data = core_contrast(rows, pixel)
            if data is None or not math.isclose(data, primary["flux_contrast"], rel_tol=1e-8, abs_tol=1e-9):
                raise ValueError(f"Could not reproduce original contrast for {obs_id}/{p['id']}")
            fitted = core_contrast(model_rows, pixel)
            residual = window(*residual_screen, primary["predicted_observed_angstrom"], p["direction"],
                              protocol["velocity_half_window_kms"])
            results[p["id"]].append({"observation_id": obs_id, "observed_core_contrast": data,
                                     "model_core_contrast": fitted,
                                     "model_to_observed_ratio": None if fitted is None else fitted/data,
                                     "residual_search": residual["outcome"]})
    summary = []
    for pred, pairs in results.items():
        ratios = [r["model_to_observed_ratio"] for r in pairs if r["model_to_observed_ratio"] is not None]
        summary.append({"id": pred, "observed_hits": len(pairs), "usable_models": len(ratios),
                        "same_sign_in_model": sum(r > 0 for r in ratios),
                        "median_model_to_observed_contrast": median(ratios) if ratios else None,
                        "residual_window_covered": sum(r["residual_search"] != "unsearchable" for r in pairs),
                        "residual_detections": sum(r["residual_search"] == "hit" for r in pairs)})
    return {"experiment": "post_holdout_template_subtraction_v1", "claim_status": "screening_only",
            "protocol_sha256": evidence["protocol_sha256"], "sample": "baseline holdout detections",
            "summary": summary, "measurements": results}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", type=Path)
    p.add_argument("holdout_evidence", type=Path)
    p.add_argument("--protocol", type=Path, default=Path("protocols/holdout-v1.json"))
    args = p.parse_args()
    try:
        protocol = args.protocol.read_bytes()
        evidence = json.loads(args.holdout_evidence.read_bytes())
        if evidence["protocol_sha256"] != digest(protocol):
            raise ValueError("Holdout evidence is from a different protocol")
        print(json.dumps(compare(evidence, json.loads(protocol), args.inputs), indent=2, allow_nan=False))
    except (ValueError, KeyError, OSError, ImportError) as exc:
        print(f"template residuals: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
