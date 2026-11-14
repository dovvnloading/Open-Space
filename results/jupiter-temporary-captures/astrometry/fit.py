"""Refit 2025 UQ146 using the three observed apparitions and propagate to Jupiter."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from initial_orbit import training_seed
from skyfit import AU_KM, DATA, DAY, EPOCH, EQUATOR_TO_ECLIPTIC, Gravity, SkyFit


def positions():
    packet = json.loads((DATA / "mpc-observations.json").read_text())[0]["ADES_DF"]
    rows = pd.DataFrame(packet).rename(
        columns={"obstime": "utc", "stn": "site", "obsid": "id"}
    )
    rows[["ra", "dec"]] = rows[["ra", "dec"]].astype(float)
    return rows[["id", "utc", "site", "ra", "dec"]]


def seed(gravity):
    orbit = json.loads((DATA / "mpc-orbit.json").read_text())[0]["mpc_orb"][0]
    if orbit["epoch_data"]["epoch"] + 2400000.5 != EPOCH:
        raise ValueError("MPC state epoch must match integration epoch.")
    cartesian = np.array(orbit["CAR"]["coefficient_values"]).reshape(2, 3)
    return (cartesian @ EQUATOR_TO_ECLIPTIC + gravity.state(0, 10, EPOCH)).ravel()


def propagate(gravity, state):
    def rhs(t, y):
        dr = gravity.planet_positions(t) - y[:3]
        force = (gravity.gm / np.linalg.norm(dr, axis=1) ** 3) @ dr
        return np.r_[y[3:], force]

    results = []
    for years in [(-43, -38), (43, 47)]:
        beg, end = [y * 365.25 for y in years]
        stop = beg if end < 0 else end
        path = solve_ivp(
            rhs,
            (0, stop),
            state,
            method="DOP853",
            rtol=2e-12,
            atol=1e-13,
            dense_output=True,
        )
        if not path.success:
            raise RuntimeError(path.message)

        def metrics(t):
            asteroid = path.sol(t).T
            jup_pos = gravity.planet_positions(t)[..., 5, :]
            jup_vel = gravity.planet_positions(t, 1)[..., 5, :]
            distance = np.linalg.norm(asteroid[..., :3] - jup_pos, axis=-1)
            rel = asteroid[..., 3:] - jup_vel
            speed = np.linalg.norm(rel, axis=-1)
            energy = (0.5 * speed**2 - gravity.gm[5] / distance) * (AU_KM / DAY) ** 2
            return distance, speed * AU_KM / DAY, energy

        time = np.linspace(beg, end, int(end - beg) + 1)
        d, v, e = metrics(time)
        roots = np.flatnonzero(np.diff(e < 0))
        if len(roots) != 2:
            results.append(
                {"negative_energy_boundaries": len(roots), "min_distance_au": d.min()}
            )
            continue
        boundaries = [
            brentq(lambda t: metrics(t)[2], time[i], time[i + 1]) for i in roots
        ]
        perijove = minimize_scalar(
            lambda t: metrics(t)[0], bounds=(beg, end), method="bounded"
        )
        t = perijove.x
        results.append(
            {
                "entry_tdb": Time(EPOCH + boundaries[0], format="jd", scale="tdb").isot,
                "exit_tdb": Time(EPOCH + boundaries[1], format="jd", scale="tdb").isot,
                "negative_energy_days": boundaries[1] - boundaries[0],
                "perijove_tdb": Time(EPOCH + t, format="jd", scale="tdb").isot,
                "distance_au": metrics(t)[0],
                "speed_km_s": metrics(t)[1],
                "energy_km2_s2": metrics(t)[2],
            }
        )
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kernel", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    gravity = Gravity(args.kernel)
    catalog_state = seed(gravity)
    rows = positions()
    fits = {}
    residuals = []
    wfst_nights = ["2025-10-27", "2025-10-28", "2025-10-29"]
    tests = ["2014", "2021", "2025", *wfst_nights]
    for name, mask in [
        ("all", np.ones(len(rows), bool)),
        *[(y, ~rows.utc.str.startswith(y)) for y in tests],
    ]:
        train = rows.loc[mask].reset_index(drop=True)
        model = SkyFit(train, gravity)
        if name == "all" or name in wfst_nights:
            start, seed_info = training_seed(model, train.utc.str[:10].to_numpy())
        else:
            start = catalog_state
            seed_info = {"source": "published MPC Cartesian state fit to all 15 positions"}
        try:
            state, _ = model.fit(start)
        except RuntimeError as exc:
            fits[name] = {"seed": seed_info, "fit_failed": str(exc)}
            continue
        tested = rows if name == "all" else rows.loc[~mask].reset_index(drop=True)
        err = SkyFit(tested, gravity).prediction(state)[0].reshape(-1, 2)
        rms = np.sqrt(np.mean(np.sum(err**2, axis=1)))
        residuals.append(
            tested.assign(fit=name, east_arcsec=err[:, 0], north_arcsec=err[:, 1])
        )
        fits[name] = {
            "seed": seed_info,
            "barycentric_state_au_day": state.tolist(),
            "test_rows": len(tested),
            "radial_rms_arcsec": rms,
            "test_max_arcsec": np.linalg.norm(err, axis=1).max(),
        }
        print(
            name,
            dict(
                (k, v) for k, v in fits[name].items() if k != "barycentric_state_au_day"
            ),
            flush=True,
        )
        if name == "all":
            fits[name]["passes"] = propagate(gravity, state)
            print(fits[name]["passes"], flush=True)
    fits["catalog_solution"] = propagate(gravity, catalog_state)
    pd.concat(residuals).to_csv(
        args.output / "fit-residuals.csv", index=False, float_format="%.7f"
    )
    (args.output / "fit-summary.json").write_text(json.dumps(fits, indent=2) + "\n")


if __name__ == "__main__":
    main()
