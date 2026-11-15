"""Measure the three-dimensional motion of 240P-B relative to 240P-A."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.linalg import norm
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import least_squares

from observations import match_components, read_mpc

HERE = Path(__file__).resolve().parent
REF = 2461000.5  # TDB Julian date at the differential state epoch
JAN1 = 2461041.5 - REF
AU_KM = 149597870.7
AU_DAY_TO_MS = AU_KM * 1000 / 86400
MU = .01720209895**2
IDENTITY = np.eye(3)


def unit(v):
    return v / norm(v, axis=-1, keepdims=True)


def tidal(r, gm):
    direction = unit(r)
    return gm / norm(r)**3 * (3 * np.outer(direction, direction) - IDENTITY)


class Motion:
    def __init__(self, path):
        table = pd.read_csv(path, dtype={"body": str})
        self.states = {
            name: CubicHermiteSpline(g.jd - REF, g[["x", "y", "z"]], g[["vx", "vy", "vz"]])
            for name, g in table.groupby("body")
        }
        self.a = self.states["90001213"]
        kwargs = dict(method="DOP853", rtol=2e-11, atol=2e-13, max_step=8, dense_output=True)
        self.before, self.after = [
            solve_ivp(self.derivative, (0., end), np.eye(9).ravel(), **kwargs).sol
            for end in (self.a.x.min(), self.a.x.max())
        ]

    def derivative(self, t, y):
        r, v = self.a(t), self.a(t, 1)
        radial, normal = unit(r), unit(np.cross(r, v))
        frame = np.column_stack((radial, np.cross(normal, radial), normal))
        h = norm(r) / 2.808
        water = .111262 * h**-2.15 * (1 + h**5.093)**-4.6142
        tensor = tidal(r, MU)
        for name, mass_ratio in [("5", 1047.3486), ("6", 3497.9)]:
            tensor += tidal(r - self.states[name](t), MU / mass_ratio)
        state = y.reshape(9, 9)
        change = np.zeros_like(state)
        change[:3] = state[3:6]
        change[3:6] = tensor @ state[:3] + water * frame @ state[6:]
        return change.ravel()

    def matrix(self, t):
        return (self.before if t < 0 else self.after)(t).reshape(9, 9)

    def design(self, obs):
        t = obs.jd.to_numpy() - REF + 69.184 / 86400
        earth = self.states["399"](t)
        tc = t - norm(self.a(t) - earth, axis=1) / 173.14463267
        line = self.a(tc) - earth
        sight = unit(line)
        east = unit(np.cross([0., 0., 1.], sight))
        north = np.cross(sight, east)
        project = np.stack([east, north], axis=1) * (206264.806247 / norm(line, axis=1))[:, None, None]
        matrices = np.array([self.matrix(time)[:3] for time in tc])
        delay = (project @ self.a(tc, 1)[:, :, None]).reshape(-1, 1)
        return (project @ matrices).reshape(-1, 9), delay

    def observable(self, r, v, epoch=JAN1):
        parent_v, parent_a = self.a(epoch, 1), self.a(epoch, 2)
        lag = (r @ parent_v) / (parent_v @ parent_v)
        dr = r @ v / norm(r)
        return {
            "separation_km": norm(r) * AU_KM,
            "distance_growth_m_s": dr * AU_DAY_TO_MS,
            "relative_speed_m_s": norm(v) * AU_DAY_TO_MS,
            "transverse_speed_m_s": norm(v - dr * unit(r)) * AU_DAY_TO_MS,
            "trailing_minutes": -lag * 1440,
            "offset_from_local_orbital_tangent_km": norm(r - lag * parent_v) * AU_KM,
            "time_shift_velocity_residual_m_s": norm(v - lag * parent_a) * AU_DAY_TO_MS,
        }


def solve(X, y):
    scale = norm(X, axis=0)
    matrix = X / scale
    initial = np.linalg.lstsq(matrix, y, rcond=None)[0]
    return least_squares(lambda b: matrix @ b - y, initial, jac=lambda b: matrix.copy(),
                         loss="soft_l1", f_scale=.45).x / scale


def measure(motion, table, check, terms):
    X, lag = motion.design(table)
    y = table[["E", "N"]].to_numpy().ravel()
    params = np.zeros(9)
    params[:terms] = solve(X[:, :terms], y)
    state = motion.matrix(JAN1) @ params
    result = motion.observable(state[:3], state[3:6])
    Xother, other_lag = motion.design(check)
    for name, res in [
        ("training", X @ params - y),
        ("other_telescope_sample", Xother @ params - check[["E", "N"]].to_numpy().ravel()),
        ("constant_time_lag", lag @ solve(lag, y) - y),
        ("constant_lag_other_sample", other_lag @ solve(lag, y) - check[["E", "N"]].to_numpy().ravel()),
    ]:
        result[name + "_median_arcsec"] = np.median(norm(res.reshape(-1, 2), axis=1))
    if terms == 6:
        evolution = {}
        for date in ("2025-10-11", "2026-01-01", "2026-03-28"):
            epoch = pd.Timestamp(date).to_julian_date() - REF
            state = motion.matrix(epoch) @ params
            details = motion.observable(state[:3], state[3:6], epoch)
            evolution[date] = {k: details[k] for k in
                               ("trailing_minutes", "offset_from_local_orbital_tangent_km")}
        result["evolution"] = evolution
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "measurements.json")
    parser.add_argument("--bootstrap", type=int, default=1000)
    args = parser.parse_args()
    data = HERE / "data"
    motion = Motion(data / "states.csv")
    positions = read_mpc(data / "astrometry.mpc")
    pairs = match_components(positions)
    mpc = pairs.groupby(["station", "night"])[["jd", "E", "N"]].median().reset_index()
    nordic = pd.read_csv(data / "not.csv")
    result = {"pairs": pairs.pairing.value_counts().to_dict(), "stations": mpc.station.nunique(),
              "station_nights": len(mpc), "nordic_nights": len(nordic)}
    for label, sample, other in [("mpc", mpc, nordic), ("nordic", nordic, mpc)]:
        result[label] = {name: measure(motion, sample, other, k) for name, k in
                         [("gravity", 6), ("radial_outgassing", 7), ("three_axis_outgassing", 9)]}
        # Resample whole sites (MPC) or whole nights (one NOT instrument).
        # All records for a selected site enter and leave together.
        X, _ = motion.design(sample)
        y = sample[["E", "N"]].to_numpy().ravel()
        if label == "mpc":
            group_rows = sample.groupby("station").indices.values()
        else:
            group_rows = [[i] for i in range(len(sample))]
        groups = [np.sort(np.r_[2 * np.asarray(g), 2 * np.asarray(g) + 1]) for g in group_rows]
        rng = np.random.default_rng(240)
        values = []
        for _ in range(args.bootstrap):
            idx = np.concatenate([groups[i] for i in rng.integers(len(groups), size=len(groups))])
            if np.linalg.matrix_rank(X[idx, :6]) != 6:
                continue
            state = motion.matrix(JAN1)[:, :6] @ solve(X[idx, :6], y[idx])
            values.append(motion.observable(state[:3], state[3:6]))
        if values:
            result[label]["gravity_bootstrap"] = pd.DataFrame(values).quantile([.025, .5, .975]).to_dict()
    b = motion.states["90001214"]
    result["jpl"] = motion.observable(b(JAN1) - motion.a(JAN1), b(JAN1, 1) - motion.a(JAN1, 1))
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
