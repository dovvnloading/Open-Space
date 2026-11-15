"""Profile sky-projected planet-chord orientations against observed Kepler-18 c/d overlaps.

Run from this directory with numpy, pandas, and scipy installed. Use --fetch to
refresh raw windows from the STScI FITS files (also requires astropy).
"""

import argparse
import itertools
from io import BytesIO
from pathlib import Path
import urllib.request

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
U1, U2 = 0.5256, 0.186
NORMALIZATION = 1 - U1 / 3 - U2 / 6
# planet order c, d, b; columns impact parameter, radius ratio, duration in hours.
GEOMETRIES = {
    "Cochran2011": np.array(
        [[0.593, 0.04549, 3.4872], [0.767, 0.05782, 3.679], [0.771, 0.01656, 2.076]]
    ),
    "DR25": np.array(
        [[0.027, 0.042627, 3.41018], [0.255, 0.051590, 3.5252], [0.741, 0.017643, 1.9811]]
    ),
}


# Rows of the endpoint sensitivity fit: 117 isolated c transits and 58 isolated d transits.
def isolated_geometry(fits, c_key, d_key):
    default = fits[(fits.u1 == 0.5256) & (fits.u2 == 0.186)].set_index("planet")
    return np.array(
        [
            default.loc[c_key, ["impact", "radius", "duration_hours"]].to_numpy(dtype=float),
            default.loc[d_key, ["impact", "radius", "duration_hours"]].to_numpy(dtype=float),
            GEOMETRIES["DR25"][2],
        ]
    )


folded = pd.read_csv(HERE / "isolated_fits.csv")
GEOMETRIES["ShortCadence"] = isolated_geometry(folded, "c", "d")


def fetch_windows():
    from astropy.io import fits

    chunks = []
    for product in pd.read_csv(HERE / "products.csv").itertuples():
        with urllib.request.urlopen(product.url, timeout=90) as response:
            with fits.open(BytesIO(response.read())) as hdus:
                data = hdus[1].data
                mask = abs(data["TIME"] - product.t0) < 0.4
                assert int(mask.sum()) == product.rows
                chunks.append(
                    pd.DataFrame(
                        {
                            "epoch": round(product.t0),
                            "file": product.file,
                            "time": data["TIME"][mask],
                            "sap": data["SAP_FLUX"][mask],
                            "pdc": data["PDCSAP_FLUX"][mask],
                            "pdc_err": data["PDCSAP_FLUX_ERR"][mask],
                            "quality": data["SAP_QUALITY"][mask],
                        }
                    )
                )
    pd.concat(chunks).to_csv(HERE / "windows.csv.gz", index=False, float_format="%.10g")


def circle_area(d, a, b):
    """Area of intersecting discs divided by pi, including containment."""
    safe = np.maximum(d, 1e-12)
    ca = np.clip((safe**2 + a * a - b * b) / (2 * safe * a), -1, 1)
    cb = np.clip((safe**2 + b * b - a * a) / (2 * safe * b), -1, 1)
    radical = np.maximum(0, (-safe + a + b) * (safe + a - b) * (safe - a + b) * (safe + a + b))
    partial = (a * a * np.arccos(ca) + b * b * np.arccos(cb) - 0.5 * np.sqrt(radical)) / np.pi
    return np.where(d >= a + b, 0.0, np.where(d <= abs(a - b), min(a, b) ** 2, partial))


def brightness(radius):
    mu = np.sqrt(np.maximum(0, 1 - radius**2))
    return (1 - U1 * (1 - mu) - U2 * (1 - mu) ** 2) / NORMALIZATION


def planet_shadow(x, impact, pr):
    radius = np.hypot(x, impact)
    return circle_area(radius, 1.0, pr) * brightness(radius) * 1e6


def build_events(geometry, column, cadence, nonlinear=None):
    windows = pd.read_csv(HERE / "windows.csv.gz")
    events = []
    for mid in pd.read_csv(HERE / "centers.csv").itertuples(index=False):
        event = windows[(windows.epoch == round(mid.c)) & (windows.quality == 0)].dropna(
            subset=[column]
        )
        if cadence == "best" and event.file.str.endswith("slc.fits").any():
            event = event[event.file.str.endswith("slc.fits")]
        else:
            event = event[event.file.str.endswith("llc.fits")]
        if len(event) < 10:
            continue
        t = event.time.to_numpy()
        flux = event[column].to_numpy()
        # Aggregate five short exposures for the fixed-template angular scan, retaining each
        # observed midpoint so that excised flags and short gaps are never interpolated.
        short = event.file.iloc[0].endswith("slc.fits")
        ref = np.median(flux[abs(t - mid.c) > 0.20])
        y = 1e6 * (flux / ref - 1)
        residual_noise = np.diff(y[abs(t - mid.c) > 0.20])
        sigma = np.median(abs(residual_noise - np.median(residual_noise))) / 0.953872
        order = np.argsort(t)
        t, y = t[order], y[order]
        if short:
            n = len(t) // 5 * 5
            times, samples = t[:n].reshape(-1, 5), y[:n].reshape(-1, 5)
            keep = np.ptp(times, axis=1) < 0.004
            grid, y = times[keep], samples[keep].mean(axis=1)
            err = sigma / np.sqrt(5)
        else:
            grid = t[:, None] + np.linspace(-0.5, 0.5, 17)[None, :] * 0.0204336 * 16 / 17
            err = sigma
        centers = np.array(mid)
        velocities = (
            2 * np.sqrt((1 + geometry[:, 1]) ** 2 - geometry[:, 0] ** 2) / (geometry[:, 2] / 24)
        )
        relative = grid - mid.c

        def model(z):
            offsets = z[:3] / 1440
            gain = z[3:6]
            shadow = sum(
                g * planet_shadow(v * (grid - t0 - delta), impact, r)
                for (impact, r, _), t0, delta, v, g in zip(
                    geometry, centers, offsets, velocities, gain
                )
            )
            dt = relative.mean(axis=1)
            return z[6] + z[7] * dt + z[8] * dt**2 - shadow.mean(axis=1)

        # Refit the additive-only null per event, allowing all timings, depth multipliers,
        # and the local quadratic to vary. A planet off-window has zero derivatives.
        z0 = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        lo = [-22] * 3 + [0.55] * 3 + [-np.inf] * 3
        hi = [22] * 3 + [1.5] * 3 + [np.inf] * 3
        fit = least_squares(
            lambda z: (y - model(z)) / err,
            z0,
            bounds=(lo, hi),
            x_scale=[3, 3, 3, 0.1, 0.1, 0.1, 100, 100, 100],
            diff_step=1e-4,
        )
        z = fit.x
        # A final white-noise estimate on five-minute (or 29-minute) samples is
        # included; broad trends are also projected out in every orientation comparison.
        residual = (y - model(z)) / err
        oot = abs(grid.mean(axis=1) - mid.c) > 0.2
        sigma_scale = max(1.0, np.median(abs(residual[oot] - np.median(residual[oot]))) * 1.4826)
        err *= sigma_scale
        jac = fit.jac
        uu, ss, _ = np.linalg.svd(jac, full_matrices=False)
        basis = uu[:, ss > ss.max() * 1e-8]
        if nonlinear is not None:
            angle, sign = nonlinear
            theta = np.radians(angle)
            (bc, rc, _), (bd, rd, _), _ = geometry

            def alternate(val):
                xc = velocities[0] * (grid - centers[0] - val[0] / 1440)
                xd = velocities[1] * (grid - centers[1] - val[1] / 1440)
                x = xd * np.cos(theta) - sign * bd * np.sin(theta)
                y2 = xd * np.sin(theta) + sign * bd * np.cos(theta)
                interior = (np.hypot(xc, bc) < 1 - rc) & (np.hypot(x, y2) < 1 - rd)
                overlap = circle_area(np.hypot(x - xc, y2 - bc), rc, rd)
                bump = (
                    overlap * brightness(np.hypot((x + xc) / 2, (y2 + bc) / 2)) * interior * 1e6
                ).mean(axis=1)
                return (y - model(val) - bump) / err

            trials = []
            for dx, dy in ((0, 0), (-3, 3), (3, -3), (-3, -3), (3, 3)):
                start = z.copy()
                start[:2] = np.clip(start[:2] + [dx, dy], -21.99, 21.99)
                trial = least_squares(
                    alternate,
                    start,
                    bounds=(lo, hi),
                    x_scale=[3, 3, 3, 0.1, 0.1, 0.1, 100, 100, 100],
                    diff_step=1e-4,
                )
                trials.append(trial)
            trial = min(trials, key=lambda f: np.dot(f.fun, f.fun))
            null = (y - model(z)) / err
            delta = np.dot(trial.fun, trial.fun) - np.dot(null, null)
            events.append((round(mid.c), delta, z.tolist(), trial.x.tolist()))
            continue
        coords = [
            v * (grid - t0 - delta / 1440) for v, t0, delta in zip(velocities, centers, z[:3])
        ]
        events.append((round(mid.c), coords, basis, (y - model(z)) / err, err))
    return events


def scan(geometry, events, angles):
    output = []
    (bc, rc, _), (bd, rd, _), _ = geometry
    for sgn, angle in itertools.product((1, -1), angles):
        theta = np.radians(angle)
        cost = power = 0.0
        for _, (xc, xd, _), basis, resid, err in events:
            # Place both central paths on the plane of the sky. bc is positive by convention.
            x = xd * np.cos(theta) - sgn * bd * np.sin(theta)
            y = xd * np.sin(theta) + sgn * bd * np.cos(theta)
            dist = np.hypot(x - xc, y - bc)
            # Require both occulting disks to lie entirely inside the stellar limb.
            # Contact cadences are cut at the integration-sample level.
            interior = (np.hypot(xc, bc) < 1 - rc) & (np.hypot(x, y) < 1 - rd)
            add = (
                circle_area(dist, rc, rd)
                * brightness(np.hypot((x + xc) / 2, (y + bc) / 2))
                * interior
                * 1e6
            ).mean(axis=1) / err
            p = add - basis @ (basis.T @ add)
            power += p @ p
            cost += p @ p - 2 * resid @ p
        output.append([sgn, angle, power, cost])
    return output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    args = p.parse_args()
    if args.fetch:
        fetch_windows()
    rows = []
    for name, geo in GEOMETRIES.items():
        for column, cadence in (("pdc", "best"), ("sap", "best"), ("pdc", "long")):
            ev = build_events(geo, column, cadence)
            for sign, angle, power, cost in scan(geo, ev, np.arange(-180, 180, 2)):
                rows.append(
                    dict(
                        geometry=name,
                        photometry=column,
                        cadence=cadence,
                        sign_d=sign,
                        angle_deg=angle,
                        expected_delta_chi2=power,
                        observed_delta_chi2=cost,
                    )
                )
    pd.DataFrame(rows).to_csv(HERE / "orientations.csv", index=False, float_format="%.5f")


if __name__ == "__main__":
    main()
