"""Replay the asteroid-pair measurements using the bundled observations and DE440s."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from encounter import encounter_grid, mpc_state, uncertainty_root
from orbits import DATA, Gravity, SkyFit


def same_time_positions(rows):
    first, second = [rows[rows.target == name] for name in ("2014 SD315", "2025 UC18")]
    matched = first.merge(second, on=["utc", "site"], suffixes=("_first", "_second"))
    angles = np.deg2rad(
        matched[["ra_first", "dec_first", "ra_second", "dec_second"]].to_numpy()
    )
    ra1, dec1, ra2, dec2 = angles.T
    hav = (
        np.sin((dec2 - dec1) / 2) ** 2
        + np.cos(dec1) * np.cos(dec2) * np.sin((ra2 - ra1) / 2) ** 2
    )
    matched["separation_deg"] = np.rad2deg(2 * np.arcsin(np.sqrt(hav)))
    return matched[["utc", "site", "id_first", "id_second", "separation_deg"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kernel", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    gravity = Gravity(args.kernel)
    rows = pd.read_csv(DATA / "astrometry.csv")
    same_time = same_time_positions(rows)
    same_time.to_csv(args.output / "same-time.csv", index=False, float_format="%.7f")
    first, second = "2014 SD315", "2025 UC18"
    names = [first, second]
    starts = {
        name: mpc_state(line, gravity)
        for name, line in zip(names, (DATA / "selected.mpc").read_text().splitlines())
    }
    # Both epoch-label files are in numbered/unnumbered MPCORB order.
    results = {}
    residuals = []
    ensembles = []
    rng = np.random.default_rng(26609)
    fitted = {}
    for name in names:
        sub = rows[rows.target == name].reset_index(drop=True)
        model = SkyFit(sub, gravity)
        state, (error, jacobian) = model.fit(starts[name])
        fitted[name] = state
        errors = error.reshape(-1, 2)
        residuals.append(
            sub.assign(east_arcsec=errors[:, 0], north_arcsec=errors[:, 1])
        )
        root = uncertainty_root(sub, jacobian)
        displacements = np.vstack((np.zeros(6), rng.normal(size=(48, 6))))
        ensembles.append(state + displacements @ root.T)
        results[name] = {
            "positions": len(sub),
            "site_nights": int((sub.site + " " + sub.utc.str[:10]).nunique()),
            "median_residual_arcsec": float(np.median(np.linalg.norm(errors, axis=1))),
            "barycentric_icrf_au_au_per_day": state.tolist(),
        }
    sub = rows[rows.target == second].reset_index(drop=True)
    results["withheld"] = {}
    for year in ["2003", "2014"]:
        excluded = sub.utc.str.startswith(year)
        train = sub[~excluded].reset_index(drop=True)
        test = sub[excluded].reset_index(drop=True)
        state, _ = SkyFit(train, gravity).fit(starts[second])
        error, _ = SkyFit(test, gravity).prediction(state)
        results["withheld"][year] = np.linalg.norm(
            error.reshape(-1, 2), axis=1
        ).tolist()
    pd.concat(residuals).to_csv(
        args.output / "residuals.csv", index=False, float_format="%.6f"
    )
    best, track = encounter_grid(gravity, *ensembles, step_days=4)
    pd.DataFrame(
        best, columns=["sampled_min_km", "lookback_years", "relative_m_s"]
    ).to_csv(args.output / "clones.csv", index=False, float_format="%.6f")
    pd.DataFrame(track, columns=["lookback_years", "nominal_separation_km"]).to_csv(
        args.output / "history.csv", index=False, float_format="%.4f"
    )
    results["nominal_4day"] = best[0].tolist()
    fine, _ = encounter_grid(gravity, ensembles[0][:1], ensembles[1][:1], step_days=2)
    results["nominal_2day"] = fine[0].tolist()
    # Shift the secondary along its orbit
    # by changing its orbital phase, preserving shape/orientation. The track is a counterfactual.
    from encounter import solar_system
    from orbits import EQUATOR_TO_ECLIPTIC

    sim = solar_system(gravity)
    pos, vel = fitted[second].reshape(2, 3) @ EQUATOR_TO_ECLIPTIC.T
    sim.add(x=pos[0], y=pos[1], z=pos[2], vx=vel[0], vy=vel[1], vz=vel[2])
    orbit = sim.particles[-1].orbit(primary=sim.particles[0])
    shifted = []
    for offset in (-5, -1, 1, 5):
        sim.add(
            primary=sim.particles[0],
            a=orbit.a,
            e=orbit.e,
            inc=orbit.inc,
            omega=orbit.omega,
            Omega=orbit.Omega,
            M=orbit.M + np.deg2rad(offset),
        )
        p = sim.particles[-1]
        shifted.append((np.array([p.xyz, p.vxyz]) @ EQUATOR_TO_ECLIPTIC).ravel())
    controls, _ = encounter_grid(
        gravity, np.tile(fitted[first], (4, 1)), np.array(shifted), step_days=4
    )
    results["mean_anomaly_controls_deg"] = {
        str(k): v.tolist() for k, v in zip((-5, -1, 1, 5), controls)
    }
    (args.output / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
