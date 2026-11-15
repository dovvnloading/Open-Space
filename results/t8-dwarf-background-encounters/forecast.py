"""Reproduce the Earth-view brown-dwarf encounter search using the frozen CSV inputs."""

import itertools
import json
from datetime import datetime, timedelta
from pathlib import Path

import erfa
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
MAS_PER_DEG = 3600000
AU_KM = 149597870.7
START, STOP = 2026.87, 2063.0
REFERENCE = 2016.0
NAME = "WISE J200050.19+362950.1"
FIVE = ["ra", "dec", "parallax", "pmra", "pmdec"]
IDS = ["2060116591502525952", "2060116595791115776"]


def unit_vectors(ra, dec):
    a, b = np.radians([ra, dec])
    sa, ca = np.sin(a), np.cos(a)
    sb, cb = np.sin(b), np.cos(b)
    return np.array([[-sa, ca, 0], [-ca * sb, -sa * sb, cb], [ca * cb, sa * cb, sb]])


def earth_position(t):
    return erfa.epv00(2451545.0, (np.asarray(t) - 2000.0) * 365.25)[1]["p"]


def grid(start, stop, step_days):
    return np.linspace(
        start, stop, int(np.ceil((stop - start) * 365.25 / step_days)) + 1
    )


def offset(row, origin):
    return (
        np.array(
            [
                ((row.ra - origin.ra + 180) % 360 - 180)
                * np.cos(np.deg2rad(origin.dec)),
                row.dec - origin.dec,
            ]
        )
        * MAS_PER_DEG
    )


def params(row, origin, epoch):
    pm = np.array([row.pmra, row.pmdec])
    xy = offset(row, origin) + pm * (REFERENCE - epoch)
    return np.r_[xy, row.parallax, pm]


def tracks(state, t, factors):
    return (
        state[..., :2, None]
        + state[..., 3:, None] * (t - REFERENCE)
        + state[..., 2, None, None] * factors
    )


def centroid_shift(separation, relative_parallax, mass):
    theta2 = 8.144 * relative_parallax * mass
    return separation * theta2 / (separation**2 + 2 * theta2)


def screen(lenses, fields):
    t = grid(START, STOP, 3)
    earth = earth_position(t)
    results = []
    paired = 0
    for lens in lenses.itertuples(index=False):
        # Reconstruct the circle used in the archive query. A finite aperture
        # limits completeness for background stars with exceptionally high PM.
        center = (2044 - lens.epoch) * np.array([lens.pmra, lens.pmdec])
        radius = 30 * np.hypot(lens.pmra, lens.pmdec) + 15000
        x = (
            ((fields.ra - lens.ra + 180) % 360 - 180)
            * np.cos(np.deg2rad(lens.dec))
            * MAS_PER_DEG
        )
        y = (fields.dec - lens.dec) * MAS_PER_DEG
        subset = fields[np.hypot(x - center[0], y - center[1]) <= radius + 5]
        factors = -unit_vectors(lens.ra, lens.dec)[:2] @ earth.T
        foreground = params(lens, lens, lens.epoch)
        for bg in subset.itertuples(index=False):
            paired += 1
            source = params(bg, lens, REFERENCE)
            missing = not np.all(np.isfinite(source[2:]))
            source[0:2] = offset(bg, lens)
            source = np.nan_to_num(source)
            source[2] = max(0.0, source[2])
            relative = source - foreground
            if relative[2] >= 0:
                continue
            xy = tracks(relative, t, factors)
            sep = np.hypot(*xy)
            idx = np.argmin(sep)
            shift = centroid_shift(sep, -relative[2], 0.05)
            if shift.max() >= 0.1:
                results.append(
                    {
                        "lens": lens.name,
                        "source_id": bg.source_id,
                        "year_at_min_sep": t[idx],
                        "sep_mas": sep[idx],
                        "max_shift_mas": shift.max(),
                        "G": bg.phot_g_mean_mag,
                        "ruwe": bg.ruwe,
                        "missing_motion_or_parallax": missing,
                    }
                )
    return pd.DataFrame(results).sort_values("max_shift_mas", ascending=False), paired


def fit_spitzer(rows, origin, train):
    xyz = rows[["x_km", "y_km", "z_km"]].values / AU_KM
    factor = -xyz @ unit_vectors(origin.ra, origin.dec)[:2].T
    t = 2000 + (rows.mjd.values - 51544.5) / 365.25
    dt = t - REFERENCE
    A = np.zeros((len(t), 2, 5))
    A[:, 0, 0] = 1
    A[:, 1, 1] = 1
    A[:, :, 2] = factor
    A[:, 0, 3] = dt
    A[:, 1, 4] = dt
    # Quoted astrometric uncertainties describe tangent-plane coordinates.
    errors = 1000 * rows[["ra_error_arcsec", "dec_error_arcsec"]].values
    positions = np.vstack([offset(r, origin) for r in rows.itertuples()])
    design = (A / errors[:, :, None])[train].reshape(-1, 5)
    data = (positions / errors)[train].ravel()
    sol = np.linalg.lstsq(design, data, rcond=None)[0]
    cov = np.linalg.inv(design.T @ design)
    return sol, cov, positions - A @ sol


def covariance(row):
    corr = np.eye(5)
    for (i, a), (j, b) in itertools.combinations(enumerate(FIVE), 2):
        corr[i, j] = corr[j, i] = getattr(row, f"{a}_{b}_corr")
    errors = np.array([getattr(row, n + "_error") for n in FIVE])
    return corr * np.outer(errors, errors)


def iso(t):
    return (
        datetime(2000, 1, 1, 12) + timedelta(days=(float(t) - 2000) * 365.25)
    ).strftime("%Y-%m-%d")


def sample_forecast(random, measured_bg, bg, lenses, t, factors, scale):
    n = len(lenses)
    background = random.multivariate_normal(measured_bg, covariance(bg) * scale**2, n)
    background[:, 2] = np.maximum(0, background[:, 2])
    values = []
    for chunk in np.array_split(background - lenses, 100):
        xy = tracks(chunk, t, factors)
        seps = np.hypot(xy[:, 0], xy[:, 1])
        idx = seps.argmin(axis=1)
        shifts = centroid_shift(seps, -chunk[:, 2, None], 0.05)
        values.extend(zip(t[idx], seps.min(axis=1), shifts.max(axis=1)))
    values = np.array(values)
    qs = np.percentile(values, [5, 50, 95], axis=0)
    return {
        "minimum_date": [iso(v) for v in qs[:, 0]],
        "sep_mas": qs[:, 1].tolist(),
        "max_shift_mas_005Msun": qs[:, 2].tolist(),
        "fraction_shift_above_02": float(np.mean(values[:, 2] > 0.2)),
    }


def forecast_models(origin, sources, fit, fit_cov):
    models = {}
    published = params(origin, origin, origin.epoch)
    old = origin.copy()
    old.ra = 300.209091
    old.dec = 36.497627
    old.parallax = 133.1
    old.pmra = 9.4
    old.pmdec = 377.9
    previous = params(old, origin, 2014.0)
    random = np.random.default_rng(20261115)
    n = 10000
    # Draw the same foreground realization for both events to retain lens covariance.
    lenses = random.multivariate_normal(fit, fit_cov, n)
    axes = unit_vectors(origin.ra, origin.dec)
    for sid, tc in zip(IDS, [2034.9, 2039.6]):
        bg = sources[sources.source_id == sid].iloc[0]
        bs = params(bg, origin, REFERENCE)
        t = grid(tc - 1.5, tc + 1.5, 0.5)
        earth = earth_position(t)
        factors = -axes[:2] @ earth.T
        # Target separation thresholds are many Einstein radii in the nominal fit.
        # Calculate the true maximum for every draw, including close approaches.
        out = {}
        measured_bg = bs.copy()
        bs[2] = max(0.0, bs[2])
        for label, fg in [
            ("census", published),
            ("2019_spitzer_overlapping", previous),
            ("epoch_refit", fit),
        ]:
            rel = bs - fg
            sep = np.hypot(*tracks(rel, t, factors))
            idx = sep.argmin()
            # The broad observing season only needs the Earth–Sun direction to a degree.
            ra_sun = -earth[idx] / np.linalg.norm(earth[idx])
            sun = np.degrees(np.arccos(np.clip(ra_sun @ axes[2], -1, 1)))
            out[label] = {
                "minimum_date": iso(t[idx]),
                "separation_mas": float(sep[idx]),
                "solar_elongation_deg": float(sun),
                "max_shift_mas_by_mass": {
                    str(m): float(centroid_shift(sep, -rel[2], m).max())
                    for m in [0.015, 0.03, 0.05, 0.07]
                },
            }
            if label == "epoch_refit":
                local = np.where((sep[1:-1] < sep[:-2]) & (sep[1:-1] < sep[2:]))[0] + 1
                out["nearby_local_minima_under_300mas"] = [
                    {"date": iso(t[j]), "separation_mas": float(sep[j])}
                    for j in local
                    if sep[j] < 300
                ]
        # A significantly negative Gaia parallax is a warning about model fit.
        # Compare the unphysical catalog propagation and fixed distant-source
        # solutions explicitly. Clipping supplies a sensitivity convention, not
        # a distance posterior. The remaining four Gaia terms stay as measured.
        for label, source_plx in [
            ("raw_gaia_parallax", measured_bg[2]),
            ("source_parallax_2mas", 2.0),
        ]:
            source = measured_bg.copy()
            source[2] = source_plx
            rel = source - fit
            sep = np.hypot(*tracks(rel, t, factors))
            idx = sep.argmin()
            out[label] = {
                "minimum_date": iso(t[idx]),
                "separation_mas": float(sep[idx]),
                "max_shift_mas_005Msun": float(
                    centroid_shift(sep, -rel[2], 0.05).max()
                ),
            }

        out["formal_5_50_95"] = sample_forecast(
            random, measured_bg, bg, lenses, t, factors, 1.0
        )
        out["triple_gaia_sigma_5_50_95"] = sample_forecast(
            random, measured_bg, bg, lenses, t, factors, 3.0
        )
        models[sid] = out
    return models


def check_ps1(rows, origin, sources):
    result = {}
    for sid, group in rows.groupby("gaia_id"):
        source = sources[sources.source_id == sid].iloc[0]
        state = params(source, origin, REFERENCE)
        state[2] = max(state[2], 0)
        usable = group[
            (group.psfQfPerfect > 0.85) & (group.raErr > 0) & (group.decErr > 0)
        ].copy()
        t = 2000 + (usable.obsTime.to_numpy() - 51544.5) / 365.25
        xy = np.vstack([offset(r, origin) for r in usable.itertuples()])
        factors = -unit_vectors(origin.ra, origin.dec)[:2] @ earth_position(t).T
        residual = xy - tracks(state, t, factors).T
        # Median by UTC observing night prevents individual same-night detection
        # clusters from dominating the summary. No Gaia-based residual clipping.
        night = (
            pd.DataFrame(
                {
                    "mjd": np.floor(usable.obsTime),
                    "x": residual[:, 0],
                    "y": residual[:, 1],
                    "year": t,
                }
            )
            .groupby("mjd")
            .median()
        )
        from scipy.stats import theilslopes

        slopes = [theilslopes(night[k], night.year, alpha=0.95) for k in ["x", "y"]]
        result[sid] = {
            "total_detections": len(group),
            "usable_detections": len(usable),
            "nights": len(night),
            "start": iso(t.min()),
            "end": iso(t.max()),
            "night_median_residual_mas": night[["x", "y"]].median().to_list(),
            "night_median_radial_residual_mas": float(
                np.median(np.hypot(night.x, night.y))
            ),
            "excess_pm_east_north_mas_yr": [
                [float(v.low_slope), float(v.slope), float(v.high_slope)]
                for v in slopes
            ],
        }
    return result


def main():
    def read(name):
        return pd.read_csv(HERE / name, dtype={"source_id": str})

    lenses, fields = read("lenses.csv"), read("gaia_fields.csv")
    origin = lenses[lenses.name == NAME].iloc[0]
    ranked, pairs = screen(lenses, fields)
    ranked.to_csv(HERE / "screen.csv", index=False, float_format="%.6f")
    rows = read("spitzer_epochs.csv")
    fit, cov, resid = fit_spitzer(rows, origin, np.ones(len(rows), dtype=bool))
    train = rows.mjd.values < 58000
    _, _, held = fit_spitzer(rows, origin, train)
    output = {
        "lenses": len(lenses),
        "unique_gaia_rows": len(fields),
        "pairs": pairs,
        "threshold_crossings": len(ranked),
        "screen_start_jyear": START,
        "screen_end_jyear": STOP,
        "epoch_refit_five_at_2016": fit.tolist(),
        "epoch_refit_errors": np.sqrt(np.diag(cov)).tolist(),
        "full_fit_axis_rmse_mas": np.sqrt(np.mean(resid**2, axis=0)).tolist(),
        "training_epochs": int(train.sum()),
        "held_out_epochs": int((~train).sum()),
        "held_out_axis_rmse_mas": np.sqrt(np.mean(held[~train] ** 2, axis=0)).tolist(),
        "panstarrs_background_checks": check_ps1(
            pd.read_csv(HERE / "ps1_detections.csv", dtype={"gaia_id": str}),
            origin,
            read("source_covariance.csv"),
        ),
        "source_forecasts": forecast_models(
            origin, read("source_covariance.csv"), fit, cov
        ),
    }
    (HERE / "forecast.json").write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
