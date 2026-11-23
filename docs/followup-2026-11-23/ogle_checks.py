"""Replay numerical audits for the OGLE follow-up packet; no new observations.

Run from any directory with numpy, scipy and astropy installed:
    python docs/followup-2026-11-23/ogle_checks.py

The fit replay is a separate command recorded in ogle.md. This script consumes
its saved result. Visibility is a geocentric solar-elongation approximation,
not a JWST ephemeris, roll/guide-star check, or approved observing schedule.
"""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import astropy
import numpy as np
import scipy
from astropy import units as u
from astropy.coordinates import GCRS, SkyCoord, get_sun
from astropy.time import Time
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "results/ogle-dark-lenses"
OUT = HERE / "ogle-evidence"
BASE_COMMIT = "972cd5169eea2b83d7407df84eb60c167b3c084e"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(exist_ok=True)
    spec = importlib.util.spec_from_file_location("ogle_pspl", SOURCE / "check/pspl.py")
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    fits = json.loads((OUT / "fits.json").read_text())
    deposited = json.loads((SOURCE / "check/out/fits.json").read_text())
    posterior = json.loads((SOURCE / "posteriors/2026-BLG-0058.json").read_text())
    manifest = json.loads((SOURCE / "sources.json").read_text())["deposited"]
    paths = ["data/ews/2026/blg-0058.dat", "data/forecast.csv",
             "posteriors/2026-BLG-0058-chains.npz", "posteriors/2026-BLG-0058.json",
             "check/out/fits.json"]
    hashes = {p: {"sha256": digest(SOURCE / p), "matches_manifest": digest(SOURCE / p) == manifest[p]}
              for p in paths}
    assert all(v["matches_manifest"] for v in hashes.values())
    assert (OUT / "baseline-replay.txt").read_bytes() == (SOURCE / "check/out/baseline.txt").read_bytes()
    fit_differences = {}
    for key in ("no_parallax", "u0+", "u0-"):
        old, new = deposited[key], fits[key]
        np.testing.assert_allclose(old["params"], new["params"], rtol=0, atol=1e-5)
        assert abs(old["chi2"] - new["chi2"]) < 1e-8
        fit_differences[key] = {
            "parameter_order": ["t0_hjd", "u0", "tE_days", "piE_north", "piE_east"],
            "parameter_differences_replay_minus_deposit": (np.array(new["params"]) - old["params"]).tolist(),
            "chi2_difference": new["chi2"] - old["chi2"],
        }
    forecast = list(csv.DictReader((SOURCE / "data/forecast.csv").open()))
    replay = list(csv.DictReader((OUT / "forecast-replay.csv").open()))
    assert forecast == replay
    calendar = {r["date"]: {"tabulated_hjd_minus_2450000": float(r["hjd_minus_2450000"]),
                              "midnight_utc_jd_minus_2450000": float(Time(r["date"]).jd - 2450000)}
                for r in forecast}
    assert all(v["tabulated_hjd_minus_2450000"] == v["midnight_utc_jd_minus_2450000"]
               for v in calendar.values())
    coord = SkyCoord(model.RA_DEC, unit=(u.hourangle, u.deg))
    times = Time("2027-01-01") + np.arange(365) * u.day
    elong = coord.transform_to(GCRS(obstime=times)).separation(get_sun(times)).deg
    visible = (elong >= 85) & (elong <= 135)
    starts = np.flatnonzero(visible & ~np.r_[False, visible[:-1]])
    ends = np.flatnonzero(visible & ~np.r_[visible[1:], False])
    windows = [{"first_daily_sample_utc": times[a].isot, "last_daily_sample_utc": times[b].isot}
               for a, b in zip(starts, ends)]
    dates = ["2027-04-01", "2027-07-01", "2027-07-15", "2027-08-15", "2027-09-01", "2028-09-01"]
    when = Time(dates)
    dn, de = model.parallax_offsets(when.jd)
    details = {}
    for key in ("u0+", "u0-"):
        p = np.array(fits[key]["params"])

        def separation_at(jd):
            n, e = model.parallax_offsets(np.atleast_1d(jd))
            tau = (jd - p[0]) / p[2] + p[3] * n + p[4] * e
            beta = -p[1] - p[3] * e + p[4] * n
            return np.sqrt(tau * tau + beta * beta)

        peak = minimize_scalar(lambda dt: float(separation_at(p[0] + dt)[0]),
                               bounds=(-10, 10), method="bounded",
                               options={"xatol": 1e-5})
        peak_jd = p[0] + peak.x
        astrometric_peak_jd = brentq(lambda dt: float(separation_at(p[0] + dt)[0] - np.sqrt(2)),
                                    50, 700) + p[0]
        sep = separation_at(when.jd)
        centroid_per_theta = sep / (sep * sep + 2)
        amp = model.magnification(p, when.jd, dn, de)
        dilution = fits[key]["fs"] * amp / (fits[key]["fs"] * amp + fits[key]["fb"])
        details[key] = {
            "fitted_t0_hjd": float(p[0]),
            "fitted_t0_calendar_approx": Time(p[0], format="jd").isot,
            "maximum_magnification_hjd": float(peak_jd),
            "maximum_magnification_calendar_approx": Time(peak_jd, format="jd").isot,
            "unblended_astrometric_maximum_hjd": float(astrometric_peak_jd),
            "unblended_astrometric_maximum_calendar_approx": Time(astrometric_peak_jd, format="jd").isot,
            "centroid_shift_for_assumed_thetaE_3p5_mas": {
                date: {"u": float(sep[i]), "unblended_mas": float(3.5 * centroid_per_theta[i]),
                       "lensed_source_fraction_of_total_light": float(dilution[i]),
                       "lensing_term_in_unresolved_centroid_mas": float(3.5 * centroid_per_theta[i] * dilution[i])}
                for i, date in enumerate(dates)},
        }
    lens_spec = importlib.util.spec_from_file_location("ogle_main_lens", SOURCE / "lens.py")
    main_lens = importlib.util.module_from_spec(lens_spec)
    lens_spec.loader.exec_module(main_lens)
    grid = Time("2026-11-01").jd - 2450000 + np.arange(800)
    offsets = main_lens.Parallax(grid, coord.ra.rad, coord.dec.rad, posterior["t0par"])
    peak_distributions = {}
    with np.load(SOURCE / "posteriors/2026-BLG-0058-chains.npz") as chains:
        for key in ("plus", "minus"):
            # candidate.py already flattened step-major, walker-minor samples
            # and saved [::5] of 40 walkers. An extra [::4] would select just
            # walkers 0 and 20. Use every stored row for this NEW calculation;
            # do not modify the previously frozen photometric forecast.
            samples = chains[key]
            tau = ((grid[None, :] - samples[:, 0, None]) / samples[:, 2, None]
                   + samples[:, 3, None] * offsets.dn + samples[:, 4, None] * offsets.de)
            beta = (samples[:, 1, None] + samples[:, 3, None] * offsets.de
                    - samples[:, 4, None] * offsets.dn)
            sep = np.hypot(tau, beta)
            crossing = (sep[:, :-1] < np.sqrt(2)) & (sep[:, 1:] >= np.sqrt(2))
            falling = (sep[:, :-1] >= np.sqrt(2)) & (sep[:, 1:] < np.sqrt(2))
            valid = crossing.any(axis=1)
            ind = crossing.argmax(axis=1)[valid]
            rows = np.flatnonzero(valid)
            fraction = (np.sqrt(2) - sep[rows, ind]) / (sep[rows, ind + 1] - sep[rows, ind])
            dates_jd = grid[ind] + fraction + 2450000
            quantiles = np.quantile(dates_jd, [0.025, 0.16, 0.5, 0.84, 0.975])
            peak_distributions[key] = {
                "method": "All saved, unweighted photometric posterior samples; first rising u=sqrt(2) crossing; daily-grid linear interpolation.",
                "sampling_limit": "The saved file contains 4800 rows from eight of the original 40 walkers (candidate.py flattened samples and saved [::5]). Samples are correlated; no extra thinning or new chain run is applied here.",
                "samples": int(len(samples)), "samples_with_crossing_in_grid": int(valid.sum()),
                "samples_below_sqrt2_at_grid_start": int((sep[:, 0] < np.sqrt(2)).sum()),
                "rising_crossing_counts_min_max": [int(crossing.sum(axis=1).min()), int(crossing.sum(axis=1).max())],
                "falling_crossing_counts_min_max": [int(falling.sum(axis=1).min()), int(falling.sum(axis=1).max())],
                "grid_start_hjd": float(grid[0] + 2450000),
                "grid_end_hjd": float(grid[-1] + 2450000),
                "grid_step_days": 1,
                "quantile_method": "numpy.quantile(method='linear')",
                "quantile_probabilities": [0.025, 0.16, 0.5, 0.84, 0.975],
                "quantiles_hjd": quantiles.tolist(),
                "calendar_quantiles_approx": Time(quantiles, format="jd").isot.tolist(),
            }
    p_mass = {}
    for key, value in posterior["solutions"].items():
        p_mass[key] = {"thetaE_for_3Msun_at_median_piE_mas": 3 * model.KAPPA * value["pie"][0],
                       "thetaE_for_3Msun_at_95pct_piE_mas": 3 * model.KAPPA * value["pie_95_upper"]}
    source_date = Time(2461215.7, format="jd").isot
    report = {
        "base_commit": BASE_COMMIT,
        "software": {"numpy": np.__version__, "scipy": scipy.__version__, "astropy": astropy.__version__},
        "input_hashes": hashes,
        "replay_checks": {"independent_fit_params_within_1e_minus_5": True, "baseline_stdout_equal": True,
                          "forecast_all_24_rows_equal": True},
        "fit_replay_differences": fit_differences,
        "forecast_calendar_check": calendar,
        "hjd_2461215p7_calendar_approx": source_date,
        "position_deg": {"ra": float(coord.ra.deg), "dec": float(coord.dec.deg)},
        "visibility_method": "Daily geocentric apparent solar separation; 85 <= elongation <= 135 deg. Not GTVT/APT.",
        "approximate_2027_jwst_windows": windows,
        "solar_elongations_deg": {},
        "astrometric_examples": details,
        "unblended_astrometric_peak_posterior": peak_distributions,
        "centroid_maximum_condition": "delta/thetaE=u/(u^2+2); derivative=(2-u^2)/(u^2+2)^2. A rising crossing of sqrt(2) changes the derivative from positive to negative.",
        "mass_scale_checks": p_mass,
        "limits": ["Calendar labels treat HJD as JD: heliocentric and time-scale corrections were not removed.",
                   "3.5 mas is an assumed angular Einstein radius, not a measurement or both-sign posterior.",
                   "The diluted lensing term omits variable flux weighting of an offset blend; joint modeling is required.",
                   "Solar-elongation windows omit JWST orbit, roll, guide stars, other constraints and telescope allocation."],
    }
    report["solar_elongations_deg"] = {
        d: float(coord.transform_to(GCRS(obstime=when[i])).separation(get_sun(when[i])).deg)
        for i, d in enumerate(dates)
    }
    (OUT / "checks.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"replay_checks": report["replay_checks"], "windows": windows,
                      "solar_elongations_deg": report["solar_elongations_deg"],
                      "mass_scale_checks": p_mass}, indent=2))


if __name__ == "__main__":
    main()
