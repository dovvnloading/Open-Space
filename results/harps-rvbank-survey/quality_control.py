"""Deterministic synthetic regression checks for the manual HARPS survey workflow."""
import os

import numpy as np
import pandas as pd

from search import Grid, INDICATORS, nightly, nuisance, search
from vet import mass


def rv_model(times, period, k, ecc, omega, tperi):
    """Independent Keplerian generator using fixed-point iteration in eccentric anomaly."""
    mean = 2 * np.pi * ((times - tperi) / period % 1)
    anomaly = mean.copy()
    for _ in range(60):
        anomaly = mean + ecc * np.sin(anomaly)
    nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(anomaly / 2),
                        np.sqrt(1 - ecc) * np.cos(anomaly / 2))
    return k * (np.cos(nu + omega) + ecc * np.cos(omega))


def regression_check(t, err, era, rv):
    design = nuisance(t, era)
    # Test the complete projected GLS calculation against ordinary least squares
    # on explicit augmented matrices, including the fibre and drift columns.
    f = 1 / np.array([10.1, 21.4, 37.19, 91.2])
    a = design / err[:, None]
    y = rv / err
    null_res = y - a @ np.linalg.lstsq(a, y, rcond=None)[0]
    reference = []
    for freq in f:
        phase = 2 * np.pi * (t - t.min()) * freq
        augmented = np.column_stack([design, np.sin(phase), np.cos(phase)]) / err[:, None]
        residual = y - augmented @ np.linalg.lstsq(augmented, y, rcond=None)[0]
        reference.append(1 - residual @ residual / (null_res @ null_res))
    grid = Grid(t, err, design, f)
    _, _, power = grid.scan(rv, return_first=True)
    np.testing.assert_allclose(power, reference, rtol=2e-9, atol=2e-11)

    # With only an era step, a systemic velocity and a drift, no power should remain.
    nuisance_only = design @ [320., .75, 70.]
    peak, _, _ = grid.scan(nuisance_only)
    np.testing.assert_allclose(peak, 0., atol=1e-12)


def main():
    if os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_EVENT_NAME") != "workflow_dispatch":
        raise SystemExit("Run this synthetic control in the manually dispatched workflow.")
    rng = np.random.default_rng(42881)
    # Distinct local nights over both fibre eras. Add jitter in observing hour.
    t = 2456200 + np.sort(rng.choice(1900, 160, replace=False)) + rng.uniform(.45, .65, 160)
    era = np.where(t < 2457169., "PRE", "POST")
    err = rng.uniform(.6, 1., len(t))
    period, k, ecc, omega, tp = 37.19, 16., .32, 1.15, t[0] - 4.
    step = 70.
    rv = (320 + step * (era == "POST") + .75 * (t - t.mean()) / 365.25
          + rv_model(t, period, k, ecc, omega, tp) + rng.normal(0, err))
    regression_check(t, err, era, rv)

    df = pd.DataFrame(dict(star="SYNTHETIC", bjd=t, rv=rv, err=err, era=era))
    for name in INDICATORS:
        df[name] = rng.normal(size=len(t))
    binned = nightly(df)
    if len(binned) != len(t) or set(binned.era) != {"PRE", "POST"}:
        raise AssertionError("Synthetic exposure dates were binned incorrectly")
    summary, peaks, injections = search("SYNTHETIC", df, 5000)
    if summary["eligibility"] != "searched" or not peaks or not injections:
        raise AssertionError("Eligible control is missing fits or injection results")
    lead = peaks[0]
    np.testing.assert_allclose(lead["period_d"], period, atol=.08, rtol=0)
    np.testing.assert_allclose(lead["k_ms"], k, atol=1., rtol=0)
    np.testing.assert_allclose(lead["ecc"], ecc, atol=.08, rtol=0)
    np.testing.assert_allclose(lead["offset_post_minus_pre_ms"], step, atol=1.5, rtol=0)
    if lead["fap"] > .01 or lead["delta_bic"] < 10:
        raise AssertionError("Injected orbit failed its blind detection thresholds")
    fit = rv_model(t, lead["period_d"], lead["k_ms"], lead["ecc"],
                   lead["omega_rad"], lead["tperi_bjd"])
    truth = rv_model(t, period, k, ecc, omega, tp)
    if np.sqrt(np.mean((fit - truth)**2)) > 1.:
        raise AssertionError("Reported orbital phase and amplitude do not reproduce the signal")
    if not any(r["period_d"] == 10 and r["k_ms"] == 100 and r["eccentricity"] == 0
               and r["recovered"] >= 6 for r in injections):
        raise AssertionError("High-amplitude 10-day injection control failed")

    # Around one solar mass, an edge-on Jupiter analogue produces about 12.5 m/s
    # over 11.86 years. These bounds distinguish MJup, MEarth and MSun units.
    edge_on = mass(11.862 * 365.25, 12.46, 0., 1.)
    low_inclination = mass(11.862 * 365.25, 12.46, 0., 1., sini=.5)
    np.testing.assert_allclose([edge_on, low_inclination], [1., 2.], rtol=.015)
    print("Synthetic controls passed for GLS power, fibre offset, eccentric orbit, "
          "injections and companion-mass units.")


if __name__ == "__main__":
    main()
