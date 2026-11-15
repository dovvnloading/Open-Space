"""Century-baseline forecast from Gaia with annual parallax and archival checks."""
import json
import warnings

import numpy as np
import pandas as pd
from astropy.time import Time
from erfa import ErfaWarning
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

from forecast import HERE, J2016, MAS_PER_DEG, covariance, event_grid, local_axes, parallax_factors


def closest(states, factors, grid):
    """Bracket the global minimum for each draw before continuous refinement."""
    dates, distances = [], []
    dt = (grid - J2016) / 365.25
    pf = factors(grid)
    for chunk in np.array_split(states, max(1, len(states) // 200)):
        x = chunk[:, 0, None] + chunk[:, 3, None]*dt + chunk[:, 2, None]*pf[:, 0]
        y = chunk[:, 1, None] + chunk[:, 4, None]*dt + chunk[:, 2, None]*pf[:, 1]
        t = grid[np.argmin(x*x + y*y, axis=1)]
        if np.any((t == grid[0]) | (t == grid[-1])):
            raise ValueError("A minimum reached the time search boundary")
        for _ in range(4):
            pos = chunk[:, :2] + chunk[:, 3:] * ((t-J2016)/365.25)[:, None] + chunk[:, 2:3]*factors(t)
            vel = chunk[:, 3:]/365.25 + chunk[:, 2:3]*factors(t, 1)
            acc = chunk[:, 2:3]*factors(t, 2)
            t -= np.sum(pos*vel, axis=1) / np.sum(vel*vel + pos*acc, axis=1)
        pos = chunk[:, :2] + chunk[:, 3:] * ((t-J2016)/365.25)[:, None] + chunk[:, 2:3]*factors(t)
        dates.extend(t)
        distances.extend(np.linalg.norm(pos, axis=1))
    return np.array(dates), np.array(distances)


def main():
    table = pd.read_csv(HERE / "gaia-2109.csv", index_col="source_id")
    lens, qso = table.loc[2124576942440058752], table.loc[2124576942438079872]
    state = event_grid(lens, qso)
    axes = local_axes(lens.ra, lens.dec)
    grid = np.arange(Time(2106., format="jyear", scale="tdb").mjd,
                     Time(2113., format="jyear", scale="tdb").mjd, 2.)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)  # ERFA flags dates past 2100
        factors = CubicSpline(grid, parallax_factors(grid, axes))
    names = ["ra", "dec", "parallax", "pmra", "pmdec"]
    cov = covariance(lens, names)
    cov[:2, :2] += covariance(qso, ["ra", "dec"])
    rng = np.random.default_rng(210905)
    samples = rng.multivariate_normal(state, cov, 20000)
    t, b = closest(np.vstack((state, samples)), factors, grid)
    years = Time(t, format="mjd", scale="tdb").jyear
    print("Stationary quasar nominal Julian year / separation (mas):", years[0], b[0])
    for name, vals in (("year", years[1:]), ("mas", b[1:])):
        print("2.5 16 50 84 97.5 percentiles", name, np.percentile(vals, [2.5,16,50,84,97.5]))
    free = state.copy()
    free[3:] -= [qso.pmra, qso.pmdec]
    free[2] -= qso.parallax
    tt, bb = closest(free[None, :], factors, grid)
    print("Taking catalog quasar PM and parallax literally (year, mas):",
          Time(tt[0], format="mjd", scale="tdb").jyear, bb[0])
    evidence = json.loads((HERE / "evidence-2109.json").read_text())
    tm = evidence["2MASS"]
    mjd = Time(float(tm["JD"]), format="jd", scale="utc").tdb.mjd
    measured = (np.array([float(tm["RAJ2000"])-lens.ra, float(tm["DEJ2000"])-lens.dec])
                * MAS_PER_DEG * [np.cos(np.radians(lens.dec)), 1])
    predicted = state[3:] * (mjd - J2016) / 365.25 + lens.parallax*parallax_factors(mjd, axes)
    print("1998 2MASS minus Gaia prediction east/north and total (mas):",
          measured-predicted, np.linalg.norm(measured-predicted))
    # A moderate radial velocity enters through changing distance. The exact
    # rectilinear perspective factor is 1/(1+v_r*dt/d); |v_r| is unmeasured.
    d_pc = 1000/lens.parallax
    for rv in (-100, 100):
        def sep(mjd):
            dt = (mjd-J2016)/365.25
            scale = 1/(1 + rv*1.022712165e-6*dt/d_pc)
            return np.linalg.norm(state[:2] + scale*(state[3:]*dt + state[2]*factors(mjd)))
        idx = np.argmin([sep(tmjd) for tmjd in grid])
        result = minimize_scalar(sep, method="bounded", bounds=(grid[idx]-3, grid[idx]+3))
        print("Radial velocity (km/s), year, mas:", rv,
              Time(result.x, format="mjd", scale="tdb").jyear, result.fun)
    for mass in (0.2, 0.35, 0.5):
        theta2 = 8.144 * mass * lens.parallax
        print("Illustrative mass (solar), brighter-image shift (mas):",
              mass, 2*theta2/(np.hypot(b[0],2*np.sqrt(theta2))+b[0]))


if __name__ == "__main__":
    main()
