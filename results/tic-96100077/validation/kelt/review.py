"""Retrospective KELT clock, calendar-season and inner-frequency controls.

Run from the KELT result directory after kelt.py has authenticated the input curves:
    OPENBLAS_NUM_THREADS=1 python review.py --raw-dir raw
If this script is elsewhere, --analysis-dir points to the directory containing kelt.py.
"""

import argparse
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.constants import c
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.time import Time
from scipy.optimize import least_squares, minimize_scalar

HERE = Path(__file__).resolve().parent
SECONDS_PER_DAY = 86400.0


def target_projection(times, ra, dec):
    """Geocenter to barycenter timing correction for a source at infinity (seconds)."""
    target = SkyCoord(ra * u.deg, dec * u.deg).cartesian.xyz.value
    xyz = get_body_barycentric("earth", times, ephemeris="builtin").xyz.to(u.m).value
    return np.einsum("ij,i->j", xyz, target) / c.to_value(u.m / u.s)


def rms(residual):
    return 1000 * np.sqrt(np.mean(residual**2))


def residuals(kelt, data, train, degree, delay, scale):
    return kelt.fit(data, train, degree, delay, scale)[0]


def annual_test(kelt, data, projection):
    """Profile the annual sine and cosine and the physical Earth-delay coefficient."""
    t = data.bjd.to_numpy() - 2457000
    angle = 2 * np.pi * t / 365.256
    annual = np.column_stack([np.sin(angle), np.cos(angle)])
    use = np.ones(len(data), dtype=bool)
    delay = kelt.gaia_delay(t)

    def fit_extra(terms, observed):
        # Extra timing delays enter with the same sign as the Gaia template.
        # Evaluate the outer ephemeris at the survey timestamps throughout this diagnostic.
        def profile(coef):
            shifted = data.assign(bjd=observed - terms @ coef / SECONDS_PER_DAY)
            return residuals(kelt, shifted, use, 2, delay, None)

        # Finite-difference offsets must exceed JD's ~40 microsecond precision.
        def jac(coef):
            vectors = np.eye(len(coef))
            step = 0.002 if terms.shape[1] == 1 else 0.25
            return np.column_stack(
                [
                    (profile(coef + step * v) - profile(coef - step * v)) / (2 * step)
                    for v in vectors
                ]
            )

        result = least_squares(
            profile, np.zeros(terms.shape[1]), jac=jac, ftol=1e-12, xtol=1e-9
        )
        if not result.success:
            raise RuntimeError(result.message)
        return result.x

    observed = data.bjd.to_numpy()
    coef = fit_extra(annual, observed)
    corrected = fit_extra(annual, observed + projection / SECONDS_PER_DAY)
    return {
        "annual_semiamplitude_s": np.linalg.norm(coef),
        "annual_correlation_with_earth_projection": np.corrcoef(
            annual @ coef, projection
        )[0, 1],
        "fitted_earth_delay_scale": fit_extra(projection[:, None], observed)[0],
        "annual_semiamplitude_after_one_correction_s": np.linalg.norm(corrected),
    }


def alias_test(kelt, data):
    """Fit four harmonics in each orientation, optimizing near daily alias frequencies."""
    t = data.bjd.to_numpy() - 2457000
    groups = [g.to_numpy() for _, g in data.groupby("orientation").groups.items()]
    y = data.flux.to_numpy()
    spacing = 1 / (5 * np.ptp(t))

    def error(freq):
        shape = kelt.basis(((t - kelt.REF) * freq) % 1)[0]
        residual = np.empty_like(y)
        for g in groups:
            coef = np.linalg.lstsq(shape[g], y[g], rcond=None)[0]
            residual[g] = y[g] - shape[g] @ coef
        return rms(residual)

    rows = []
    for offset in [0, -1, 1, -1.0027379, 1.0027379, -2.0054758, 2.0054758]:
        center = 1 / kelt.PERIOD + offset
        grid = np.arange(center - 0.0015, center + 0.0015, spacing)
        best = grid[np.argmin([error(f) for f in grid])]
        result = minimize_scalar(
            error,
            bounds=(best - spacing, best + spacing),
            method="bounded",
            options={"xatol": 1e-11},
        )
        rows.append(
            {
                "frequency_offset_per_day": offset,
                "best_frequency_per_day": result.x,
                "rms_ppt": result.fun,
            }
        )
    return rows


def diagnostics(kelt, data, t, year, night, degree):
    delay = kelt.gaia_delay(t)
    observations = np.ones(len(data), dtype=bool)
    unique_nights, index = np.unique(night, return_inverse=True)
    folds = index % 5
    null = np.empty(len(data))
    fixed = null.copy()
    for f in range(5):
        test = folds == f
        null[test] = residuals(kelt, data, ~test, degree, delay, 0)[test]
        fixed[test] = residuals(kelt, data, ~test, degree, delay, kelt.SCALE)[test]
    _, free, _ = kelt.fit(data, observations, degree, delay, None)
    blocks = (night - night.min()) // 21
    low, high = kelt.paired_interval(null, fixed, blocks)
    year_low, year_high = kelt.paired_interval(null, fixed, year)
    gains = null**2 - fixed**2
    summary = {
        "degree": degree,
        "heldout_nights": len(unique_nights),
        "zero_night_rms_ppt": rms(null),
        "fixed_night_rms_ppt": rms(fixed),
        "gain_night_mse_ppt2": rms(null) ** 2 - rms(fixed) ** 2,
        "block21_gain_ci025_ppt2": low,
        "block21_gain_ci975_ppt2": high,
        "year_gain_ci025_ppt2": year_low,
        "year_gain_ci975_ppt2": year_high,
        "fitted_full_scale": free[-1],
    }
    rows = []
    for y in np.unique(year):
        test = year == y
        null = residuals(kelt, data, ~test, degree, delay, 0)[test]
        fixed = residuals(kelt, data, ~test, degree, delay, kelt.SCALE)[test]
        _, free, _ = kelt.fit(data, ~test, degree, delay, None)
        rows.append(
            {
                "degree": degree,
                "withheld_year": y,
                "interleaved_gain_this_year_ppt2": np.mean(gains[test]) * 1e6,
                "n": test.sum(),
                "zero_year_rms_ppt": rms(null),
                "fixed_year_rms_ppt": rms(fixed),
                "gain_year_mse_ppt2": rms(null) ** 2 - rms(fixed) ** 2,
                "fitted_scale_excluding_year": free[-1],
            }
        )
    return summary, rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--analysis-dir", type=Path, default=HERE)
    p.add_argument("--raw-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, default=HERE)
    args = p.parse_args()
    sys.path.insert(0, str(args.analysis_dir.resolve()))
    import kelt
    from measure import DEC, RA

    inputs = kelt.load(args.raw_dir, False)
    tables = {"clock": [], "aliases": [], "calibration": [], "omitted-seasons": []}
    for reduction, data in inputs.groupby("reduction"):
        data = data.reset_index(drop=True)
        observed = data.bjd.to_numpy()
        time = Time(observed, format="jd", scale="tdb")
        projection = target_projection(time, RA, DEC)
        # Preserve the observing-night folds for every comparison and clock convention.
        nights = np.floor(observed).astype(int)
        years = np.asarray(time.ymdhms.year)
        row = annual_test(kelt, data, projection)
        tables["clock"].append(dict(reduction=reduction, **row))
        tables["aliases"].extend(
            dict(reduction=reduction, **r) for r in alias_test(kelt, data)
        )
        for name, offset in [
            ("archive", 0),
            ("plus_earth", projection / SECONDS_PER_DAY),
        ]:
            shifted = data.assign(bjd=observed + offset)
            t = shifted.bjd.to_numpy() - 2457000
            for degree in [1, 2]:
                result, seasons = diagnostics(kelt, shifted, t, years, nights, degree)
                tables["calibration"].append(
                    dict(reduction=reduction, clock=name, **result)
                )
                tables["omitted-seasons"].extend(
                    dict(reduction=reduction, clock=name, **s) for s in seasons
                )
    args.output.mkdir(exist_ok=True, parents=True)
    for name, rows in tables.items():
        pd.DataFrame(rows).to_csv(
            args.output / f"review-{name}.csv", index=False, float_format="%.6f"
        )


if __name__ == "__main__":
    main()
