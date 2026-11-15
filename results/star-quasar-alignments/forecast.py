"""Recompute the 2028 geocentric encounter and independent imaging checks."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import get_body_barycentric
from astropy.time import Time
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
MAS_PER_DEG = 3_600_000.0
J2016 = Time(2016.0, format="jyear", scale="tdb").mjd


def local_axes(ra, dec):
    a, d = np.radians([ra, dec])
    return np.array([[-np.sin(a), np.cos(a), 0],
                     [-np.cos(a) * np.sin(d), -np.sin(a) * np.sin(d), np.cos(d)]])


def parallax_factors(mjd, axes):
    earth = get_body_barycentric("earth", Time(mjd, format="mjd", scale="tdb"))
    return -(axes @ earth.xyz.to_value("AU")).T


def covariance(row, names):
    sigmas = row[[name + "_error" for name in names]].to_numpy(float)
    corr = np.eye(len(names))
    for i, a in enumerate(names):
        for j, b in enumerate(names[:i]):
            corr[i, j] = corr[j, i] = row[b + "_" + a + "_corr"]
    return corr * np.outer(sigmas, sigmas)


def event_grid(lens, source):
    # RA errors and Gaia pmra are already multiplied by cos(dec).
    scale = MAS_PER_DEG * np.cos(np.radians(lens.dec))
    dx = (lens.ra - source.ra) * scale
    dy = (lens.dec - source.dec) * MAS_PER_DEG
    return np.array([dx, dy, lens.parallax, lens.pmra, lens.pmdec])


def solve(events, factors):
    """Vectorized minimization on the single close approach around March 2028."""
    t = np.full(len(events), Time("2028-03-20", scale="tdb").mjd)
    for _ in range(5):
        pos = (events[:, :2] + events[:, 3:] * ((t - J2016) / 365.25)[:, None]
               + events[:, 2:3] * factors(t))
        vel = events[:, 3:] / 365.25 + events[:, 2:3] * factors(t, 1)
        acc = events[:, 2:3] * factors(t, 2)
        t -= np.sum(pos * vel, axis=1) / np.sum(vel * vel + pos * acc, axis=1)
    pos = (events[:, :2] + events[:, 3:] * ((t - J2016) / 365.25)[:, None]
           + events[:, 2:3] * factors(t))
    return t, np.linalg.norm(pos, axis=1)


def independent_imaging(lens, source, axes, generator):
    det = pd.read_csv(HERE / "ps1.csv")
    det = det[(det.psfFlux > 5 * det.psfFluxErr) & (det.psfFluxErr > 0)
              & (det.psfQfPerfect > 0.85) & (det.infoFlag2 == 0)].copy()
    mjd = det.obsTime.to_numpy()
    year = (mjd - J2016) / 365.25
    det["x"] = (det.ra - source.ra) * MAS_PER_DEG * np.cos(np.radians(lens.dec))
    det["y"] = (det.dec - source.dec) * MAS_PER_DEG
    positions = det[["x", "y"]].to_numpy()
    initial = event_grid(lens, source)
    motion = parallax_factors(mjd, axes) * lens.parallax
    lens_xy = initial[:2] + initial[3:] * year[:, None] + motion
    # Windows assign the two well-separated PS1 images. All centroids within
    # either 1 arcsecond window enter the fits, with no clipping on residuals.
    det["night"] = np.floor(mjd).astype(int)
    det["dt"] = year
    q = det[np.linalg.norm(positions, axis=1) < 1000].groupby("night")[["x", "y"]].median()
    det[["x", "y"]] -= motion
    s = det[np.linalg.norm(positions - lens_xy, axis=1) < 1000]
    s = s.groupby("night")[["dt", "x", "y"]].median()

    def ps_solution(star, qso):
        design = np.column_stack((np.ones(len(star)), star.dt))
        pos, pm = np.linalg.lstsq(design, star[["x", "y"]], rcond=None)[0]
        return [*(pos - qso[["x", "y"]].median()), lens.parallax, *pm]

    state = np.array(ps_solution(s, q))
    # Per-object night bootstrap conditional on Gaia parallax. This diagnostic
    # does not include common frame or spatially varying systematics.
    def samples(frame):
        return frame.iloc[generator.integers(len(frame), size=len(frame))]
    draws = np.array([ps_solution(samples(s), samples(q)) for _ in range(1500)])
    # Use a global daily grid for the less precise ground-based trajectories.
    # Some draws have broad minima which defeat local Newton iteration.
    grid = np.arange(Time("2023-01-01", scale="tdb").mjd,
                     Time("2034-01-01", scale="tdb").mjd)
    pf = parallax_factors(grid, axes)
    dt = (grid - J2016)/365.25
    t, u = [], []
    for v in np.vstack((state, draws)):
        xy = v[:2] + dt[:, None]*v[3:] + v[2]*pf
        rad = np.linalg.norm(xy, axis=1)
        index = rad.argmin()
        if index in (0, len(rad)-1):
            raise ValueError("PS1 minimum on search boundary")
        t.append(grid[index])
        u.append(rad[index])
    t, u = np.asarray(t), np.asarray(u)
    print("PS1 nights, foreground/background:", len(s), len(q))
    print("PS1 proper motion (mas/yr):", np.round(state[3:], 2))
    print("PS1 central date and separation (mas):", utc(t[0]), round(u[0], 2))
    print("PS1 bootstrap date and separation percentiles (16, 50, 84):",
          [utc(x)[:10] for x in np.percentile(t[1:], [16, 50, 84])],
          np.percentile(u[1:], [16, 50, 84]).round(2))
    print("PS1 quasar offset from Gaia (mas, east/north):", q.median().to_numpy().round(2))


def utc(mjd):
    return Time(mjd, format="mjd", scale="tdb").utc.isot


def main():
    table = pd.read_csv(HERE / "gaia.csv", index_col="source_id")
    lens, source = table.loc[1018659143306715392], table.loc[1018659143306715008]
    axes = local_axes(lens.ra, lens.dec)
    day = np.arange(Time("2026-01-01", scale="tdb").mjd,
                    Time("2030-01-01", scale="tdb").mjd + 1)
    factors = CubicSpline(day, parallax_factors(day, axes))
    state = event_grid(lens, source)
    t0, d0 = solve(state[None, :], factors)
    # Independent scalar minimization checks the stationarity solver.
    def separation(mjd):
        return np.linalg.norm(state[:2] + state[3:] * (mjd - J2016) / 365.25
                              + state[2] * factors(mjd))
    check = minimize_scalar(separation, bounds=(t0[0]-100, t0[0]+100), method="bounded")
    assert abs(check.x - t0[0]) < 0.01 and abs(check.fun - d0[0]) < 0.001
    cov = covariance(lens, ["ra", "dec", "parallax", "pmra", "pmdec"])
    cov[:2, :2] += covariance(source, ["ra", "dec"])
    rng = np.random.default_rng(11032028)
    states = rng.multivariate_normal(state, cov, size=50000)
    time, dist = solve(states, factors)
    print("Gaia nominal closest approach (UTC, mas):", utc(t0[0]), round(d0[0], 3))
    print("Gaia 2.5, 16, 50, 84, 97.5 percentiles UTC:",
          [utc(x)[:10] for x in np.percentile(time, [2.5, 16, 50, 84, 97.5])])
    print("Gaia separation percentiles (mas):", np.percentile(dist, [2.5, 16, 50, 84, 97.5]).round(2))
    model = json.loads((HERE / "literature-values.json").read_text())
    for label, mass in (("H", model["wd"]["MassH"]), ("He", model["wd"]["MassHe"]),
                        ("2022 assumed", model["kluter"]["mass"])):
        theta2 = 8.144 * mass * lens.parallax
        u = d0[0] / np.sqrt(theta2)
        major = 2 * theta2 / (np.hypot(d0[0], 2*np.sqrt(theta2)) + d0[0])
        amplification = (u*u + 2) / (u * np.sqrt(u*u + 4))
        fr = 10**(-0.4*(source.phot_g_mean_mag - lens.phot_g_mean_mag))
        shift = d0[0] / (u*u + 2)
        blend = fr*amplification*(d0[0]+shift)/(1+fr*amplification) - fr*d0[0]/(1+fr)
        print(label, "mass, Einstein radius, major-image shift, blended G shift (mas),",
              "quasar brightening (mmag):",
              np.round([mass, np.sqrt(theta2), major, blend, 2500*np.log10(amplification)], 4))
    independent_imaging(lens, source, axes, rng)
    lines = pd.read_csv(HERE / "sdss-lines.csv")
    lines["snr"] = lines.LINEAREA / lines.LINEAREA_ERR
    print(lines[(lines.LINEAREA_ERR > 0) & (lines.snr > 3)][["LINENAME","LINEZ","snr"]].to_string(index=False))


if __name__ == "__main__":
    main()
