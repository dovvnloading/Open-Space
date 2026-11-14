"""Predict an independent recovery of 2025 UC18 near its next opposition."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time
from scipy.integrate import solve_ivp

from encounter import uncertainty_root
from orbits import ARCSEC, DATA, EPOCH, Gravity, LIGHT_AU_DAY, SkyFit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kernel", type=Path)
    parser.add_argument("replay_directory", type=Path)
    args = parser.parse_args()
    gravity = Gravity(args.kernel)
    summary = json.loads((args.replay_directory / "summary.json").read_text())
    rows = pd.read_csv(DATA / "astrometry.csv")
    selected = rows[rows.target == "2025 UC18"].reset_index(drop=True)
    state = np.array(summary["2025 UC18"]["barycentric_icrf_au_au_per_day"])
    _, jacobian = SkyFit(selected, gravity).prediction(state)
    root = uncertainty_root(selected, jacobian)
    times = Time(
        pd.date_range("2026-11-15", "2027-12-31", freq="D").to_numpy(), scale="utc"
    )
    jd = times.tdb.jd
    earth = gravity.state(0, 3, jd)[0] + gravity.state(3, 399, jd)[0]
    sun = gravity.state(0, 10, jd)[0] - earth
    solution = solve_ivp(
        gravity.acceleration,
        (0, jd[-1] - EPOCH),
        np.r_[state, np.eye(6).ravel()],
        method="DOP853",
        atol=5e-13,
        rtol=5e-12,
        dense_output=True,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    emission = jd - EPOCH
    for _ in range(4):
        results = solution.sol(emission).T
        vector = results[:, :3] - earth
        distance = np.linalg.norm(vector, axis=1)
        emission = jd - EPOCH - distance / LIGHT_AU_DAY
    direction = vector / distance[:, None]
    elongation = np.rad2deg(
        np.arccos(
            np.clip(
                np.sum(direction * sun, axis=1) / np.linalg.norm(sun, axis=1), -1, 1
            )
        )
    )
    center = np.argmax(elongation)
    records = []
    for k in [center - 7, center, center + 7]:
        lon = np.arctan2(direction[k, 1], direction[k, 0])
        east = np.array([-np.sin(lon), np.cos(lon), 0])
        axes = np.vstack((east, np.cross(direction[k], east)))
        deriv = results[k, 6:].reshape(6, 6)[:3]
        velocity = results[k, 3:6]
        deriv = deriv - np.outer(velocity, direction[k] @ deriv) / (
            LIGHT_AU_DAY + direction[k] @ velocity
        )
        sigmas = np.linalg.svd(
            axes @ deriv @ root * ARCSEC / distance[k], compute_uv=False
        )
        records.append(
            [
                times[k].utc.isot,
                np.rad2deg(lon) % 360,
                np.rad2deg(np.arcsin(direction[k, 2])),
                distance[k],
                elongation[k],
                *sigmas,
            ]
        )
    pd.DataFrame(
        records,
        columns=[
            "geocentric_utc",
            "icrf_ra_deg",
            "icrf_dec_deg",
            "range_au",
            "solar_elongation_deg",
            "formal_major_arcsec",
            "formal_minor_arcsec",
        ],
    ).to_csv(args.replay_directory / "prediction.csv", index=False, float_format="%.8f")


if __name__ == "__main__":
    main()
