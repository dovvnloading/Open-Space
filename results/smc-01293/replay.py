"""Refit OGLE-IV only, project frozen waveforms onto external observations."""
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
BANDS = ("g", "r", "i", "z")
TARGET = "0122m730.2381"
EPOCH = 2457000.0


def read_csv(name):
    with gzip.open(DATA / name, "rt") as stream:
        return list(csv.DictReader(stream))


def clean(data):
    time, mag, error = data.T
    return data[np.isfinite(time + mag + error) & (error > 0) &
                (error < max(2.5 * np.median(error), 0.017))]


def season(time):
    return np.floor((time - 2455100) / 365.25).astype(int)


def harmonics(time, period):
    angle = 2 * np.pi * (time - EPOCH) / period
    return np.column_stack([f(k * angle) for k in range(1, 9) for f in (np.cos, np.sin)])


def regress(design, values, errors, nights, noise_floor):
    sigma = np.hypot(errors, noise_floor)
    weights = 1 / sigma
    values = values - np.median(values)
    for step in range(4):
        coeff = np.linalg.lstsq(design * weights[:, None], values * weights, rcond=None)[0]
        residual = values - design @ coeff
        if step < 3:
            weights = np.sqrt(np.minimum(1, 2.5 * sigma / np.maximum(abs(residual), 1e-8))) / sigma
    z = abs(residual / sigma)
    loss = np.sum(np.where(z < 2.5, z**2, 5*z - 6.25))
    bread = np.linalg.pinv(design.T @ (weights[:, None]**2 * design))
    influence = bread @ (design.T * (weights**2 * residual))
    groups, inverse = np.unique(nights, return_inverse=True)
    by_night = np.zeros((len(groups), design.shape[1]))
    np.add.at(by_night, inverse, influence.T)
    return coeff, np.sqrt((by_night**2).sum(axis=0)), loss


def indicators(labels):
    return [labels == value for value in np.unique(labels)]


def template_fit(data, periods, extra=None):
    time, mag, err = data.T
    columns = indicators(season(time))
    for p in periods:
        columns.extend(harmonics(time, p).T)
    if extra is not None:
        columns.extend(extra.T)
    return regress(np.column_stack(columns), mag, err, time.astype(int), 0.004)


def fit_periods(train):
    initial = np.array([0.6703844, 0.529915537])
    def objective(delta):
        return template_fit(train, initial + np.asarray(delta) * 1e-6)[2]
    optimized = minimize(objective, [0, 0], method="Nelder-Mead",
                         options={"initial_simplex": [[0, 0], [0.1, 0], [0, 0.1]],
                                  "xatol": 0.0007, "fatol": 0.003})
    if not optimized.success:
        raise RuntimeError(optimized.message)
    periods = initial + optimized.x * 1e-6
    coeff, _, loss = template_fit(train, periods)
    return periods, coeff[-32:].reshape((2, 16)), loss


def waveforms(time, periods, coeffs):
    return [harmonics(time, p) @ c for p, c in zip(periods, coeffs)]


def project(data, periods, coeffs, labels=None, filter_errors=True, floor=0.005):
    if filter_errors:
        data = clean(data)
    time, mag, errors = data.T
    if labels is None:
        labels = season(time)
    design = np.column_stack([*indicators(labels), *waveforms(time, periods, coeffs)])
    coeff, se, _ = regress(design, mag, errors, time.astype(int), floor)
    return {"points": len(time), "nights": len(np.unique(time.astype(int))),
            "primary_scale": round(float(coeff[-2]), 4), "primary_se": round(float(se[-2]), 4),
            "secondary_scale": round(float(coeff[-1]), 4), "secondary_se": round(float(se[-1]), 4)}


def get_photometry():
    datasets = {}
    for r in read_csv("photometry.csv.gz"):
        key = r["dataset"]
        datasets.setdefault(key, []).append([float(r["hjd"]), float(r["mag"]), float(r["err"])])
    return {key: np.asarray(rows) for key, rows in datasets.items()}


def heliocentric_mjd(mjd, seconds):
    jd = mjd + 2400000.5 + seconds / 172800
    days = jd - 2451545
    lon = np.deg2rad((280.460 + 0.9856474 * days) % 360)
    anomaly = np.deg2rad((357.528 + 0.9856003 * days) % 360)
    lon += np.deg2rad(1.915 * np.sin(anomaly) + 0.020 * np.sin(2 * anomaly))
    obl = np.deg2rad(23.439 - 0.0000004 * days)
    radius = 1.00014 - 0.01671*np.cos(anomaly) - 0.00014*np.cos(2*anomaly)
    ra, dec = np.deg2rad([11.964375, -73.1199722])
    project_sun = (np.cos(lon)*np.cos(ra)*np.cos(dec) +
                   np.sin(lon)*np.cos(obl)*np.sin(ra)*np.cos(dec) +
                   np.sin(lon)*np.sin(obl)*np.sin(dec))
    return jd - radius * project_sun * 499.005 / 86400


def decam(periods, coeffs):
    rows = read_csv("decam.csv.gz")
    centers = {}
    for obj in {r["objid"] for r in rows}:
        for band in BANDS:
            samples = [r for r in rows if r["objid"] == obj and r["filter"] == band]
            if samples:
                centers[obj, band] = np.median([[float(r[k]) for k in ("mag", "ra", "dec")] for r in samples], axis=0)
    measurements = {}
    for r in rows:
        offset = np.array([float(r[k]) for k in ("mag", "ra", "dec")]) - centers[r["objid"], r["filter"]]
        measurements.setdefault(r["expnum"], {})[r["objid"]] = (r, offset)
    grouped = {b: [] for b in BANDS}
    for objs in measurements.values():
        if TARGET not in objs or len(objs) < 4:
            continue
        r, delta = objs[TARGET]
        ref = np.median([v for obj, (_, v) in objs.items() if obj != TARGET], axis=0)
        hjd = heliocentric_mjd(float(r["mjd"]), float(r["seconds"]))
        m, dx, dy = delta - ref
        grouped[r["filter"]].append([hjd, m, float(r["err"]), dx*3600*np.cos(np.deg2rad(-73.1199722)), dy*3600])
    out = {}
    for band, pts in grouped.items():
        arr = np.array(pts)
        out[band] = project(arr[:, :3], periods, coeffs, np.zeros(len(arr)), False, floor=0.01)
    pts = np.vstack([v for v in grouped.values()])
    bands = np.concatenate([[b]*len(v) for b,v in grouped.items()])
    out["all_bands"] = project(pts[:, :3], periods, coeffs, bands, False, floor=0.01)
    # Coefficients in arcsec/mag. The two templates are fitted simultaneously.
    for i, axis in ((3,"east"),(4,"north")):
        series = np.column_stack([pts[:, 0], pts[:, i], np.full(len(pts), 0.04)])
        out[axis] = project(series, periods, coeffs, bands, False, floor=0)
    return out


def analyze():
    series = get_photometry()
    train = clean(series["OGLE-IV-I"])
    periods, coeffs, loss = fit_periods(train)
    p2 = periods[1]
    angle = 2*np.pi*(train[:, 0] - EPOCH) / (2*p2)
    odd = np.column_stack([f(k*angle) for k in (1, 3, 5, 7, 9) for f in (np.cos, np.sin)])
    odd_delta = loss - template_fit(train, periods, odd)[2]
    grid = np.linspace(EPOCH, EPOCH + 20, 50000)
    depth = [float(np.ptp(v)) for v in waveforms(grid, periods, coeffs)]
    alias_loss = {}
    f0 = 1 / periods[1]
    for offset in (-1.00274, -1., 1., 1.00274, 2.00548):
        def objective(delta):
            q = 1 / (f0 + offset + delta)
            return template_fit(train, (periods[0], q))[2]
        optimum = minimize_scalar(objective, bounds=(-0.00013, 0.00013), method="bounded",
                                  options={"xatol": 1e-8})
        alias_loss[str(offset)] = round(float(optimum.fun - loss), 2)
    ii = clean(series["OGLE-II-I"])
    iii = clean(series["OGLE-III-I"])
    common = np.isin(iii[:, 0], ii[:, 0])
    unique_iii = iii[~common]
    out = {"ogle_shared_hjd": int(common.sum()),
           "ogle_iii_i_excluding_ii_epochs": project(unique_iii, periods, coeffs, filter_errors=False),
           "period_days": periods.tolist(), "peak_to_peak_mag": depth,
           "double_period_loss_improvement_10_terms": float(odd_delta),
           "daily_alias_excess_robust_loss": alias_loss,
           "survey_projections": {key: project(s, periods, coeffs) for key, s in series.items()},
           "decam": decam(periods, coeffs)}
    return out, periods, coeffs


if __name__ == "__main__":
    print(json.dumps(analyze()[0], indent=2))
