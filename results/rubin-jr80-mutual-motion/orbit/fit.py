"""Search periods and component quadrants of an astrometric binary orbit."""
import argparse
import csv
import itertools
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks

BASE = Path(__file__).resolve().parent
PARAMETERS = "freq_per_day eccentricity phase semimajor_at40au_arcsec node cos_i omega".split()
rows = list(csv.DictReader((BASE / "observations.csv").open()))


def column(key):
    return np.array([float(row[key]) for row in rows])


t, east, north, ra, dec, dist, sigma, block = map(
    column, "mjd_utc east_arcsec north_arcsec ra_deg dec_deg distance_AU sigma_arcsec block".split()
)
t -= dist * 499.004783836 / 86400  # differential light travel time
reference = [row["id"] for row in rows].index("HST")
t -= t[reference]
sky = np.column_stack((east, north))
block = block.astype(int)


def component_signs(independent_rubin):
    """Allow each separated historical observing block to reverse its component labels."""
    groups = block.copy()
    if independent_rubin:
        groups[-3:] = [5, 6, 7]
    free = [group for group in np.unique(groups) if group != 4]
    options = []
    for choices in itertools.product((-1, 1), repeat=len(free)):
        signs_by_group = dict(zip(free, choices))
        options.append([signs_by_group.get(group, 1) for group in groups])
    return np.asarray(options)


signs = component_signs(False)


def tangent_axes(ra_deg, dec_deg):
    ra, dec = np.deg2rad(ra_deg), np.deg2rad(dec_deg)
    east = np.stack((-np.sin(ra), np.cos(ra), np.zeros_like(ra)), axis=1)
    north = np.stack((-np.cos(ra)*np.sin(dec), -np.sin(ra)*np.sin(dec), np.cos(dec)), axis=1)
    return east, north


u, v = tangent_axes(ra, dec)
basis = np.array([u[reference], v[reference], np.cross(u[reference], v[reference])])
scale = 40 / dist


def ell(mean_anomaly, eccentricity):
    """Kepler-plane coordinates in units of the semimajor axis."""
    anomaly = mean_anomaly + eccentricity * np.sin(mean_anomaly)
    for _ in range(7):
        anomaly -= (anomaly-eccentricity*np.sin(anomaly)-mean_anomaly)/(1-eccentricity*np.cos(anomaly))
    return np.stack((np.cos(anomaly)-eccentricity,
                     np.sqrt(1-eccentricity**2)*np.sin(anomaly)), axis=-1)


def initial(mask, frequencies):
    """At each trial phase and eccentricity, solve four linear sky-plane coefficients."""
    error = sigma[mask]
    measurements = signs[:, mask, None] * (sky[mask] / error[:, None])
    stacked = np.moveaxis(measurements, 0, -1).reshape(len(error), -1)
    power = np.sum(measurements**2, axis=(1, 2))
    profile = np.full((len(frequencies), 5), np.nan)
    profile[:, 0] = np.inf
    angles = 2*np.pi*frequencies[:, None]*t[None, mask]
    for e in np.arange(0, .901, .1):
        for phase in np.arange(0, 2*np.pi, .4 if e > 0 else 7):
            design = ell(angles+phase, e) * (scale[mask]/error)[None, :, None]
            transposed = design.transpose(0, 2, 1)
            cross = transposed @ stacked
            coefficients = np.linalg.solve(transposed @ design, cross)
            improvement = (cross*coefficients).reshape(len(frequencies), 2, 2, -1).sum(axis=(1, 2))
            cost = power - improvement
            quadrant = np.argmin(cost, axis=1)
            minimum = cost[np.arange(len(frequencies)), quadrant]
            better = minimum < profile[:, 0]
            profile[better] = np.column_stack((minimum[better], frequencies[better],
                                               np.full(better.sum(), e), np.full(better.sum(), phase),
                                               quadrant[better]))
    return profile


def predict(parameters, time=None, east_axis=None, north_axis=None, angular_scale=None):
    frequency, e, phase, a, node, cos_i, omega = parameters
    first = np.array([np.cos(node), np.sin(node), 0.])
    second = np.array([-np.sin(node)*cos_i, np.cos(node)*cos_i, np.sqrt(max(0, 1-cos_i**2))])
    periapsis = (first*np.cos(omega) + second*np.sin(omega)) @ basis
    orthogonal = (-first*np.sin(omega) + second*np.cos(omega)) @ basis
    offsets = a * ell(2*np.pi*frequency*(t if time is None else time)+phase, e) @ np.stack((periapsis, orthogonal))
    east = np.sum(offsets*(u if east_axis is None else east_axis), axis=1)
    north = np.sum(offsets*(v if north_axis is None else north_axis), axis=1)
    return (scale if angular_scale is None else angular_scale)[:, None]*np.column_stack((east, north))


def refine(seed, mask):
    _, f, e, phase, quadrant = seed
    quadrant = int(quadrant)
    signed = signs[quadrant, :, None] * sky
    design = ell(2*np.pi*f*t+phase, e) * scale[:, None]
    coefficients = np.linalg.lstsq(design[mask]/sigma[mask, None],
                                     signed[mask]/sigma[mask, None], rcond=None)[0].T
    left, axes, right = np.linalg.svd(coefficients)
    if np.linalg.det(left) < 0:
        left[:, 1] *= -1
        axes[1] *= -1
    if np.linalg.det(right) < 0:
        right[1, :] *= -1
        axes[1] *= -1
    node = np.arctan2(left[1, 0], left[0, 0])
    omega = np.arctan2(right[1, 0], right[0, 0])

    def residual(parameters):
        return ((predict(parameters)[mask]-signed[mask])/sigma[mask, None]).ravel()

    alternatives = []
    for flip in (0, np.pi):  # reflected orbit orientations
        start = [f, e, phase, axes[0], node+flip, np.clip(axes[1]/axes[0], -.99, .99), omega+flip]
        fit = least_squares(residual, start, max_nfev=120, ftol=1e-8,
                            bounds=([f-.0007, 0, -30, .05, -30, -.9999, -30],
                                    [f+.0007, .94, 30, 10, 30, .9999, 30]))
        alternatives.append((2*fit.cost, fit.x, quadrant))
    return min(alternatives, key=lambda result: result[0])


def perform(mask, output):
    frequencies = np.arange(1/1100, 1/14, .000012)
    chunks = np.array_split(frequencies, max(1, len(signs)//18))
    grid = np.concatenate([initial(mask, chunk) for chunk in chunks])
    local = find_peaks(-grid[:, 0], distance=3)[0]
    seeds = sorted(local, key=lambda i: grid[i, 0])[:200]
    results = sorted((refine(grid[i], mask) for i in seeds), key=lambda result: result[0])
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["chi2","period_days",*PARAMETERS,"signindex","heldout_axial_error_arcsec"])
        for cost, parameters, quadrant in results:
            observed = sky[~mask]
            predicted = predict(parameters)[~mask]
            error = np.minimum(np.linalg.norm(predicted-observed, axis=1),
                               np.linalg.norm(predicted+observed, axis=1))
            writer.writerow([cost, 1/parameters[0], *parameters, quadrant, " ".join(f"{v:.5f}" for v in error)])


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--omit", default="all", help="epoch ID, or 'no_rubin' to omit all three Rubin nights")
    cli.add_argument("--independent-rubin", action="store_true", help="separate quadrant choices at every Rubin epoch")
    args = cli.parse_args()
    signs = component_signs(args.independent_rubin)
    mask = np.ones(len(rows), dtype=bool)
    if args.omit == "no_rubin":
        mask[-3:] = False
    elif args.omit in ("2020", "2021"):
        season = np.array([r["id"].startswith("CFHT_") and int(float(r["mjd_utc"])) // 100 == (590 if args.omit == "2020" else 593) for r in rows])
        if args.omit == "2021":
            season |= np.array([r["id"] == "CFHT_59403" for r in rows])
        mask &= ~season
    elif args.omit != "all":
        ids = [row["id"] for row in rows]
        if args.omit not in ids:
            cli.error("unknown epoch ID")
        mask[ids.index(args.omit)] = False
    qualifier = "quadrants_" if args.independent_rubin else ""
    perform(mask, BASE / f"fit_{qualifier}{args.omit}.csv")
