"""Propagate the full, conditional A1/A2 fits into the proposed observing windows."""
import argparse
from importlib import util
import json
import os
from pathlib import Path
import tempfile

import numpy as np


HERE = Path(__file__).resolve().parent


def compute(build, orbits, settings):
    spec = util.spec_from_file_location("orbit_replay", orbits / "replay.py")
    engine = util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    # The forward adapter expects stations.txt beside its input pack.
    engine.HERE = HERE
    os.environ["FO_ACCEL"] = "1"
    states = json.loads((HERE / "force_states.json").read_text())
    gravity = np.array(json.loads((orbits / "fits/VE138.json").read_text())["full"]["state"])
    sky = []
    for r in settings:
        sky.append({"utc": r["utc"], "site": r["site"], "mode": "CCD", "catalog": "Gaia3",
                    "ra_deg": str(r["gravity_ra_deg"]), "dec_deg": str(r["gravity_dec_deg"])})

    with tempfile.TemporaryDirectory() as work:
        model = engine.Forward(build, Path(work), sky, "2026 VV15",
                               str(states["epoch_jd_tdb"]), gravity)
        try:
            origin = model(np.r_[gravity, 0.0, 0.0])
            if np.max(np.abs(origin)) > 0.002:
                raise RuntimeError("Gravity reference does not reproduce the frozen pointings")
            offsets = {key: model(np.array(state)) - origin
                       for key, state in states["model_states"].items()}
            full = np.array(states["model_states"]["A12"])
            matrix = np.asarray(states["A12_site_utc_date_cluster_covariance"])
            sigma = np.sqrt(np.diag(matrix))
            # Scale to dimensionless correlations before diagonalization.
            weights, axes = np.linalg.eigh(matrix / np.outer(sigma, sigma))
            projected = []
            for size, axis in zip(weights, axes.T):
                perturb = np.sqrt(max(0, size)) * axis * sigma
                projected.append((model(full + perturb) - model(full - perturb)) / 2)
        finally:
            model.close()
    forecasts = []
    for i, row in enumerate(settings):
        item = dict(row)
        directions = np.array(projected)[:, i, :]
        variances, vectors = np.linalg.eigh(directions.T @ directions)
        fits = {}
        for key, xy in offsets.items():
            east, north = xy[i]
            dec = row["gravity_dec_deg"]
            fits[key] = {"ra_deg": round(row["gravity_ra_deg"] + east / 3600 / np.cos(np.radians(dec)), 9),
                         "dec_deg": round(dec + north / 3600, 9),
                         "offset_east_north_arcsec": [round(east, 4), round(north, 4)],
                         "offset_tangent_plane_arcsec": round(np.hypot(east, north), 4)}
        item["force_models"] = fits
        item["A12_cluster_sigma_major_minor_arcsec"] = np.sqrt(variances[::-1]).round(4).tolist()
        item["A12_major_axis_pa_deg_east_of_north"] = round(np.degrees(np.arctan2(vectors[0, 1], vectors[1, 1])) % 180, 3)
        forecasts.append(item)
    return forecasts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("--orbits", type=Path, default=HERE.parent / "orbits")
    parser.add_argument("--output", type=Path, help="Write a regenerated result instead of checking it")
    args = parser.parse_args()
    path = HERE / "force_forecasts.json"
    reference = json.loads(path.read_text())
    settings = [{k: v for k, v in row.items() if k not in
                 ("force_models", "A12_cluster_sigma_major_minor_arcsec", "A12_major_axis_pa_deg_east_of_north")}
                for row in reference]
    result = compute(args.build.resolve(), args.orbits.resolve(), settings)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    else:
        for expected, actual in zip(reference, result, strict=True):
            for key, fit in expected["force_models"].items():
                if not np.allclose(actual["force_models"][key]["offset_east_north_arcsec"],
                                   fit["offset_east_north_arcsec"], atol=0.04):
                    raise RuntimeError(f"Force projection mismatch at {actual['utc']}")
            if not np.allclose(actual["A12_cluster_sigma_major_minor_arcsec"],
                               expected["A12_cluster_sigma_major_minor_arcsec"], atol=0.08):
                raise RuntimeError(f"Covariance projection mismatch at {actual['utc']}")
            print(f"{actual['utc']} {actual['site']}: A1+2 vs gravity "
                  f"{actual['force_models']['A12']['offset_tangent_plane_arcsec']:.2f} arcsec")


if __name__ == "__main__":
    main()
