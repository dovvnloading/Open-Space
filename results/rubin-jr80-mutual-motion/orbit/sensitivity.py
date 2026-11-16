"""Refit each leading alias under alternative PSF measurements and error budgets."""
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

import fit

BASE = Path(__file__).resolve().parent


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def nightly_vectors(case):
    measured = fit.sky.copy()
    if case not in ("diff_only", "empirical_psf"):
        return measured
    output = BASE.parent / "output"
    visits = {r["id"]: r for r in read(output / "astrometry.csv")}
    if case == "diff_only":
        stamps = [r for r in read(output/"fits.csv") if r["kind"]=="Difference" and r["variant"]=="nominal"]
        axes = ["east_arcsec", "north_arcsec"]
    else:
        stamps = [r for r in read(output/"empirical_fits.csv") if r["mode"]=="target_all"]
        axes = ["east", "north"]
    nights = {}
    for row in stamps:
        if row["id"] in visits:
            obs = visits[row["id"]]
            vec = np.array([float(row[key]) for key in axes])
            vec *= np.sign(vec @ np.array([float(obs["east"]), float(obs["north"])]))
            nights.setdefault(obs["utc_date"], []).append(vec)
    for j, epoch in enumerate(fit.rows):
        if epoch["id"] in nights:
            measured[j] = np.mean(nights[epoch["id"]], axis=0)
    return measured


def physics(x):
    period = 1 / x[0]
    axis_km = x[3] * 40 * 149597870.7 / 206264.806247
    mass = 4*np.pi**2 * (axis_km*1000)**3 / (6.67430e-11 * (period*86400)**2)
    return np.array([period, x[1], axis_km, mass])


def checks(case, starting_models):
    sky = nightly_vectors(case)
    sigma = fit.sigma.copy()
    mask = np.ones(len(sigma), dtype=bool)
    if case == "inflate_ground":
        sigma[sigma > .01] *= 1.75
    if case == "block_weights":
        for group in np.unique(fit.block):
            member = fit.block == group
            sigma[member] *= np.sqrt(member.sum())
    if case == "no_rubin":
        mask[-3:] = False

    results = []
    for record in starting_models:
        start = np.array([float(record[key]) for key in fit.PARAMETERS])
        f = start[0]
        signs = fit.signs[int(record["signindex"]), :, None]

        def residual(x):
            return ((fit.predict(x)-signs*sky)[mask]/sigma[mask,None]).ravel()

        solution = least_squares(residual, start,
                                 bounds=([f-.00035, .001, -20, .07, -30, -.999, -30],
                                         [f+.00035, .96, 30, 10, 30, .999, 30]))
        cov = np.linalg.pinv(solution.jac.T @ solution.jac)
        dx = 1e-7 * np.eye(7)
        derived_jac = np.array([(physics(solution.x+d)-physics(solution.x-d))/(2e-7) for d in dx]).T
        if np.linalg.matrix_rank(solution.jac) == 7:
            errors = np.sqrt(np.maximum(0, np.diag(derived_jac @ cov @ derived_jac.T)))
        else:
            errors = np.full(4, np.nan)
        results.append([case, 2*solution.cost, *physics(solution.x), *errors])
    return sorted(results, key=lambda row: row[1])[:5]


if __name__ == "__main__":
    initial = read(BASE/"fit_all.csv")[:12]
    with (BASE/"sensitivity.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow("case chi2 period_days eccentricity semimajor_km total_mass_kg se_P se_e se_a_km se_M_kg".split())
        for choice in ["nominal", "inflate_ground", "block_weights", "diff_only", "no_rubin", "empirical_psf"]:
            writer.writerows(checks(choice, initial))
