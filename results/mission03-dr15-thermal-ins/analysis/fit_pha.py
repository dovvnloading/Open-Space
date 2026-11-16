"""Fit a PPS source/background PHA using its ARF and RMF and joint Poisson counts."""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq
from scipy.special import xlogy
from astropy.io import fits

# Morrison & McCammon (1983), ApJ 270, 119, table 2. sigma in cm2 per H atom.
BREAKS = np.array(
    [
        0.03,
        0.1,
        0.284,
        0.4,
        0.532,
        0.707,
        0.867,
        1.303,
        1.840,
        2.471,
        3.210,
        4.038,
        7.111,
        8.331,
    ]
)
C0 = np.array(
    [
        17.3,
        34.6,
        78.1,
        71.4,
        95.5,
        308.9,
        120.6,
        141.3,
        202.7,
        342.7,
        352.2,
        433.9,
        629.0,
        701.2,
    ]
)
C1 = np.array(
    [
        608.1,
        267.9,
        18.8,
        66.8,
        145.8,
        -380.6,
        169.3,
        146.8,
        104.7,
        18.7,
        18.7,
        -2.4,
        30.9,
        25.2,
    ]
)
C2 = np.array(
    [
        -2150.0,
        -476.1,
        4.3,
        -51.4,
        -61.1,
        294.0,
        -47.7,
        -31.5,
        -17.0,
        0.0,
        0.0,
        0.75,
        0.0,
        0.0,
    ]
)


def sigma(e):
    i = np.maximum(0, np.searchsorted(BREAKS, e, side="right") - 1)
    return 1e-24 * (C0[i] + e * C1[i] + e * e * C2[i]) / e**3


def photon_shape(e, shape, param):
    if shape == "blackbody":
        return e**2 / np.expm1(np.minimum(e / param, 700))
    if shape == "brems_continuum":
        return np.exp(-e / param) / e
    return e ** (-param)


def background_expectation(n, b, alpha, s):
    # Profile nuisance expectation beta in background aperture; n~Pois(s+a*beta), b~Pois(beta).
    linear = (1 + alpha) * s - alpha * (n + b)
    root = np.sqrt(linear**2 + 4 * alpha * (1 + alpha) * b * s)
    beta = (root - linear) / (2 * alpha * (1 + alpha))
    # Rationalize the positive root where subtraction loses precision.
    np.divide(2 * b * s, root + linear, out=beta, where=linear > 0)
    return np.maximum(beta, 1.0e-200)


def profile_stat(n, b, alpha, s):
    beta = background_expectation(n, b, alpha, s)
    mu = np.maximum(s + alpha * beta, 1.0e-200)
    return 2 * np.sum(mu - n + xlogy(n, n / mu) + beta - b + xlogy(b, b / beta))


def read(path, minimum):
    with fits.open(path) as f:
        h = f["SPECTRUM"]
        meta = h.header.copy()
        chan = h.data["CHANNEL"].astype(int)
        n = h.data["COUNTS"].astype(float)
        good = h.data["QUALITY"] == 0
    with fits.open(path.parent / meta["BACKFILE"]) as f:
        h = f["SPECTRUM"]
        b = h.data["COUNTS"].astype(float)
        if not np.array_equal(chan, h.data["CHANNEL"]):
            raise ValueError("background channel grid mismatch")
        alpha = (
            meta["EXPOSURE"]
            * meta["BACKSCAL"]
            / (h.header["EXPOSURE"] * h.header["BACKSCAL"])
        )
    p = path.parent / meta["RESPFILE"]
    if not p.exists():
        p = p.with_name(p.name + ".gz")
    with fits.open(p) as f:
        rmf = f["MATRIX"].data
        eb = f["EBOUNDS"].data
        true_lo, true_hi = (
            np.asarray(rmf["ENERG_LO"], float),
            np.asarray(rmf["ENERG_HI"], float),
        )
        e = (true_hi + true_lo) / 2
        if not np.array_equal(eb["CHANNEL"], chan):
            raise ValueError("PHA/RMF channels differ")
        eb_mid = (eb["E_MIN"] + eb["E_MAX"]) / 2
        bounds = np.unique(
            np.round(
                np.r_[
                    np.arange(minimum, 1.001, 0.05),
                    np.arange(1.1, 2.001, 0.1),
                    2.5,
                    3,
                    4,
                    5,
                ],
                6,
            )
        )
        bins = np.searchsorted(bounds, eb_mid, side="right") - 1
        nb = len(bounds) - 1
        legit = good & (bins >= 0) & (bins < nb)
        sn = np.bincount(bins[legit], weights=n[legit], minlength=nb)
        bn = np.bincount(bins[legit], weights=b[legit], minlength=nb)
        m = np.zeros((nb, len(e)))
        for i, row in enumerate(rmf):
            k = int(row["N_GRP"])
            if not k:
                continue
            starts = np.atleast_1d(row["F_CHAN"])[:k].astype(int)
            lengths = np.atleast_1d(row["N_CHAN"])[:k].astype(int)
            c = np.concatenate([np.arange(s, s + z) for s, z in zip(starts, lengths)])
            channel_idx = c - chan[0]
            mask = legit[channel_idx]
            vals = np.asarray(row["MATRIX"])[: len(c)]
            m[:, i] = np.bincount(
                bins[channel_idx[mask]], weights=vals[mask], minlength=nb
            )
    with fits.open(path.parent / meta["ANCRFILE"]) as f:
        a = f["SPECRESP"].data
        if not (
            np.allclose(true_lo, a["ENERG_LO"]) and np.allclose(true_hi, a["ENERG_HI"])
        ):
            raise ValueError("ARF and RMF energy grids differ")
        m *= a["SPECRESP"] * (true_hi - true_lo) * meta["EXPOSURE"]
    # Drop bins with no good channels. Ignore out-of-calibration true energies below 0.03 keV.
    keep = np.bincount(bins[legit], minlength=nb) > 0
    te = (e >= 0.03) & (e <= 10.0)
    return (
        e[te],
        m[keep][:, te],
        sn[keep],
        bn[keep],
        float(alpha),
        bounds[:-1][keep],
        bounds[1:][keep],
        true_lo[te],
        true_hi[te],
    )


def profile_shape(shape, m, nh_scale, n, b, alpha, limits, grid_points):
    """Minimize over source counts and bounded absorbing column at a fixed shape."""
    low, high = np.exp(limits[0])
    high = min(high, max(low, float(n.sum())))

    def at_column(nh):
        counts = m @ (np.exp(-nh * nh_scale) * shape)
        weights = counts / counts.sum()

        def score(total):
            source = total * weights
            mu = source + alpha * background_expectation(n, b, alpha, source)
            return 1.0 - np.sum(weights * (n / mu))

        # At fixed shape the profiled Poisson likelihood is concave in source counts.
        if score(low) >= 0:
            total = low
        elif score(high) <= 0:
            total = high
        else:
            total = brentq(score, low, high)
        return profile_stat(n, b, alpha, total * weights)

    if grid_points < 3:
        raise ValueError("Column profile requires at least three grid points")
    nh_low, nh_high = limits[1]
    if nh_low == nh_high:
        return at_column(nh_low)
    # Concentrate trials near small columns, and refine every sampled local minimum.
    grid = nh_low + (nh_high - nh_low) * np.linspace(0, 1, grid_points) ** 2
    values = [at_column(nh) for nh in grid]
    candidates = list(values)
    for i, value in enumerate(values):
        left, right = max(i - 1, 0), min(i + 1, len(grid) - 1)
        if value <= values[left] and value <= values[right]:
            optimum = minimize_scalar(
                at_column,
                bounds=(grid[left], grid[right]),
                method="bounded",
                options={"xatol": 1.0e-9},
            )
            if not optimum.success or not np.isfinite(optimum.fun):
                raise RuntimeError("Absorbing-column profile did not converge")
            candidates.append(optimum.fun)
    return min(candidates)


def run(path, minimum, fixed_nh=None, *, profile_grid_points=9):
    energy, m, n, b, alpha, lo, hi, _, _ = read(path, minimum)
    nh_scale = sigma(energy) * 1.0e21
    result = {
        "file": path.name,
        "energy_edges_keV": [[float(l), float(h)] for l, h in zip(lo, hi)],
        "source_counts": n.astype(int).tolist(),
        "background_counts": b.astype(int).tolist(),
        "background_scale": alpha,
        "models": {},
    }
    for model, shape_bounds in [
        ("blackbody", [0.02, 0.8]),
        ("powerlaw", [0.5, 8.0]),
        ("brems_continuum", [0.025, 10.0]),
    ]:
        limits = [
            (-8, 16),
            (0, 10) if fixed_nh is None else (fixed_nh, fixed_nh),
            tuple(np.log(shape_bounds)),
        ]

        def stat(x):
            total, nh, theta = np.exp(x[0]), x[1], np.exp(x[2])
            shape = photon_shape(energy, model, theta)
            unit = m @ (np.exp(-nh * nh_scale) * shape)
            counts = total * unit / max(np.sum(unit), 1e-200)
            return profile_stat(n, b, alpha, counts)

        starts = np.geomspace(*shape_bounds, 5)[1:-1]
        fits_ = [
            minimize(
                stat,
                [np.log(max(1.0, n.sum() - alpha * b.sum())), nh, np.log(st)],
                bounds=limits,
                method="L-BFGS-B",
            )
            for st in starts
            for nh in [0.0, 1.0]
        ]
        converged = [fit for fit in fits_ if fit.success and np.isfinite(fit.fun)]
        if not converged:
            raise RuntimeError(f"No converged {model} fit")
        best = min(converged, key=lambda fit: fit.fun)
        total, nh, theta = np.exp(best.x[0]), best.x[1], np.exp(best.x[2])
        norm = total / np.sum(
            m @ (np.exp(-nh * nh_scale) * photon_shape(energy, model, theta))
        )

        def delta(theta_log):
            shape = photon_shape(energy, model, np.exp(theta_log))
            return (
                profile_shape(
                    shape, m, nh_scale, n, b, alpha, limits, profile_grid_points
                )
                - best.fun
            )

        if abs(delta(best.x[2])) > 1.0e-5:
            raise RuntimeError(f"{model} minimum disagrees with the profiled deviance")
        ci = []
        endpoint_deltas = []
        for end in limits[-1]:
            endpoint = None
            achieved = None
            if abs(end - best.x[2]) > 1.0e-6 and delta(end) > 2.706:
                root = brentq(lambda z: delta(z) - 2.706, *sorted((best.x[2], end)))
                achieved = float(delta(root))
                if abs(achieved - 2.706) > 1.0e-5:
                    raise RuntimeError(f"{model} profile endpoint failed verification")
                endpoint = float(np.exp(root))
            ci.append(endpoint)
            endpoint_deltas.append(achieved)
        egrid = np.linspace(0.2, 2, 4001)
        energy_flux = 1.602176634e-9 * np.trapezoid(
            norm
            * photon_shape(egrid, model, theta)
            * np.exp(-nh * 1e21 * sigma(egrid))
            * egrid,
            x=egrid,
        )
        result["models"][model] = {
            "deviance": float(best.fun),
            "normalization": norm,
            "NH_1e21": nh,
            "shape": theta,
            "shape_profile90": ci,
            "shape_profile90_delta_deviance": endpoint_deltas,
            "absorbed_flux_0.2_2": float(energy_flux),
            "optimizer_converged": bool(best.success),
        }
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("spectrum", type=Path)
    p.add_argument("--minimum", type=float, choices=(0.2, 0.3), default=0.2)
    p.add_argument(
        "--fixed-nh", type=float, help="Freeze column in units of 1e21 cm^-2"
    )
    p.add_argument(
        "--profile-grid-points",
        type=int,
        default=9,
        help="Number of absorbing-column search points before local refinement",
    )
    a = p.parse_args()
    result = run(
        a.spectrum, a.minimum, a.fixed_nh, profile_grid_points=a.profile_grid_points
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
