"""Replay the catalog search and encounter-distance tests from the frozen tables."""
import csv
import hashlib
import json
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
PAR = ["ra", "dec", "parallax", "pmra", "pmdec", "radial_velocity"]
K = 1.022712165      # pc/Myr per km/s
C = 299792.458       # km/s
G = 0.00449850215    # pc^3 / (solar mass Myr^2)


def load(name):
    with (ROOT / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def vizier_row(name):
    # VizieR TSV has a header, units, ruler and data, interleaved with metadata.
    lines = [s for s in (ROOT / name).read_text().splitlines()
             if s.strip() and not s.startswith("#")]
    assert len(lines) == 4
    return dict(zip(lines[0].split("\t"), lines[-1].split("\t"), strict=True))


def parameters(rows):
    return np.array([[float(row[key]) for key in PAR] for row in rows])


def state(parameters, light_time=False):
    ra, dec, parallax, mu_east, mu_north, rv = parameters.T
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    cd, sd, ca, sa = np.cos(dec), np.sin(dec), np.cos(ra), np.sin(ra)
    radial = np.column_stack((cd * ca, cd * sa, sd))
    east = np.column_stack((-sa, ca, np.zeros_like(ra)))
    north = np.column_stack((-sd * ca, -sd * sa, cd))
    distance = 1000.0 / parallax
    position = radial * distance[:, None]
    tangential = 4.74047046 / parallax[:, None] * (
        mu_east[:, None] * east + mu_north[:, None] * north)
    if light_time:
        tangential *= (1 + rv / C)[:, None]
    velocity = K * (rv[:, None] * radial + tangential)
    if light_time:
        position += velocity * (distance / (C * K))[:, None]
    return position, velocity


def closest(dp, dv):
    v2 = np.sum(dv * dv, axis=-1)
    t = -np.sum(dp * dv, axis=-1) / v2
    return np.linalg.norm(dp + np.expand_dims(t, -1) * dv, axis=-1), t, np.sqrt(v2) / K


def point(par, light_time=False):
    p, v = state(par, light_time)
    return closest(p[1] - p[0], v[1] - v[0])


def search():
    rows = load("hosts.csv")
    p, v = state(parameters(rows))
    accepted = 0
    ranked = []
    for i, star in enumerate(rows[:-1]):
        d, t, speed = closest(p[i + 1:] - p[i], v[i + 1:] - v[i])
        valid = (speed > 5) & (abs(t) < 2)
        accepted += int(valid.sum())
        for k in np.flatnonzero(valid & (d < 1)):
            ranked.append((d[k], t[k], speed[k], star["hostname"], rows[i + 1 + k]["hostname"]))
    ranked.sort()
    assert [ranked[0][3], ranked[0][4]] == ["HD 34445", "HD 37124"]
    return dict(hosts=len(rows), total_pairs=len(rows) * (len(rows) - 1) // 2,
                pass_time_and_speed=accepted, inside_1pc=len(ranked),
                inside_0p1pc=sum(1 for x in ranked if x[0] < 0.1),
                closest_five=ranked[:5])


def uncertain_state(row, rng, n):
    sigma = np.array([float(row[key + "_error"]) for key in PAR])
    cov = np.diag(sigma * sigma)
    for i, j in combinations(range(5), 2):
        correlation = float(row[PAR[i] + "_" + PAR[j] + "_corr"])
        cov[i, j] = cov[j, i] = correlation * sigma[i] * sigma[j]
    noise = rng.multivariate_normal(np.zeros(6), cov, n)
    # Gaia ra_error is an angular offset in alpha* = alpha cos(delta).
    noise[:, 0] /= 3_600_000 * np.cos(np.deg2rad(float(row["dec"])))
    noise[:, 1] /= 3_600_000
    return state(parameters([row]) + noise)


def galactic_model(par):
    # ICRS to Galactic Cartesian (x toward Galactic centre). The Sun then has x<0.
    rotation = np.array([[-.0548755604, -.8734370902, -.4838350155],
                         [.4941094279, -.4448296300, .7469822445],
                         [-.8676661490, -.1980763734, .4559837762]])
    pos, vel = state(par)
    p = pos @ rotation.T + [-8200, 0, 20.8]
    v = vel @ rotation.T + K * np.array([11.1, 245.0, 7.25])

    def derivative(t, state_vector):
        r, speed = state_vector.reshape(2, 2, 3)
        x, y, z = r.T
        # Miyamoto–Nagai disk and two spherical Plummer components, pc and Msun.
        beta = np.hypot(z, 250.0)
        a_beta = 3500 + beta
        f = -G * 7.91e10 / (x*x + y*y + a_beta*a_beta)**1.5
        acc = np.column_stack((x, y, a_beta * z / beta)) * f[:, None]
        for mass, scale in [(1.4e10, 350), (69.8e10, 24000)]:
            acc -= G * mass * r / (np.sum(r*r, axis=1) + scale*scale)[:, None]**1.5
        return np.r_[speed.ravel(), acc.ravel()]

    solution = solve_ivp(derivative, [0, .22], np.r_[p.ravel(), v.ravel()],
                         rtol=1e-12, atol=1e-12, dense_output=True, method="DOP853")
    assert solution.success
    def distance(t):
        r = solution.sol(t)[:6].reshape(2, 3)
        return np.linalg.norm(r[1] - r[0])
    opt = minimize_scalar(distance, bounds=(.19, .22), method="bounded",
                          options={"xatol": 1e-14})
    assert opt.success
    return float(opt.fun), float(opt.x)


def follow_up():
    pair = sorted(load("gaia_pair.csv"), key=lambda x: x["source_id"])
    nominal = parameters(pair)
    keck = nominal.copy()
    keck[:, 5] = [float(r["mean_rv_kms"]) for r in load("keck.csv")]
    long_pm, hip = keck.copy(), keck.copy()
    hip_only_pm = keck.copy()
    for i, hip_id in enumerate([24681, 26381]):
        hg = vizier_row(f"hgca_{hip_id}.tsv")
        hp = vizier_row(f"hip_{hip_id}.tsv")
        assert hg["Gaia"].strip() == pair[i]["source_id"]
        long_pm[i, 3:5] = [float(hg[k]) for k in ("pmRAhg", "pmDEhg")]
        hip[i, :5] = [float(hp[k]) for k in ("RArad", "DErad", "Plx", "pmRA", "pmDE")]
        hip_only_pm[i, 3:5] = hip[i, 3:5]
    dr2 = parameters(sorted(load("gaia_dr2.csv"), key=lambda x: x["source_id"]))
    q = [0.05, 0.5, 0.95]
    rng = np.random.default_rng(37124)
    n = 200000
    (pa, va), (pb, vb) = [uncertain_state(row, rng, n) for row in pair]
    d, t, v = closest(pb - pa, vb - va)
    result = {"nominal_pc_myr_kms": point(nominal),
              "current_distance_pc": float(np.linalg.norm(np.diff(state(nominal)[0], axis=0))),
              "gaia_mc_draws": n, "quantiles": q, "gaia_distance_pc": np.quantile(d, q).tolist(),
              "gaia_time_myr": np.quantile(t, q).tolist(), "gaia_speed_kms": np.quantile(v, q).tolist(),
              "gaia_fraction_inside_0p1pc": float(np.mean(d < .1)),
              "gaia_dr2_epoch2015p5": point(dr2),
              "gaia_positions_keck_rv": point(keck), "hg_long_pm_keck_rv": point(long_pm),
              "hip_only_pm_gaia_parallax_keck_rv": point(hip_only_pm),
              "hip_astrometry_keck_rv_epoch1991p25": point(hip),
              "galactic_tide_pc_myr": galactic_model(nominal),
              "light_time_pc_myr_kms": point(nominal, light_time=True)}
    for label, index, shift, base in [("rv_offsets_0p5_kms", 5, .5, nominal),
                                      ("parallax_offsets_0p03_mas", 2, .03, nominal)]:
        grid = []
        for offsets in product([-shift, 0, shift], repeat=2):
            trial = base.copy()
            trial[:, index] += offsets
            dist, time, _ = point(trial)
            grid.append([*offsets, dist, time])
        result[label] = grid
    return result


if __name__ == "__main__":
    for filename, expected in json.loads((ROOT / "provenance.json").read_text())["frozen_files"].items():
        if hashlib.sha256((ROOT / filename).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Input digest mismatch: {filename}")
    print(json.dumps({"search": search(), "follow_up": follow_up()}, indent=2))
