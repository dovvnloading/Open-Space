"""Inspect the eclipse-depth leads, using nights as resampling units."""

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

STARS = ["07730", "06703", "09564", "04349", "07125", "00863"]
EDGES = (0, 2457000, 2458300, 2459700, 2500000)
ERAS = ("2010–14", "2014–18", "2018–20", "2022–24")


def phase(t, epoch, period):
    return ((t - epoch) / period + 0.5) % 1 - 0.5


def depth_measure(t, p, m, center, rng):
    distance = abs((p - center + 0.5) % 1 - 0.5)
    core = distance < 0.025
    reference = (abs(p) > 0.16) & (abs(p) < 0.34)
    selected = core | reference
    # Count each local observing night once per window. Use HJD integer boundaries
    # (noon UT, daylight in Chile), with the same resampling for the two windows.
    nights, inv = np.unique(np.floor(t[selected]), return_inverse=True)
    if len(nights) == 0:
        return None
    groups = np.full((2, len(nights)), np.nan)
    for i, mask in enumerate([core[selected], reference[selected]]):
        for j in np.unique(inv[mask]):
            groups[i, j] = np.median(m[selected][mask & (inv == j)])
    counts = np.isfinite(groups).sum(axis=1)
    if min(counts) < 4:
        return None
    med = np.nanmedian(groups, axis=1)
    # Both sets can be sparse in V. Resample all nights jointly, dropping
    # replicate draws with an empty group.
    samples = rng.integers(0, len(nights), size=(1800, len(nights)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        bs = np.nanmedian(groups[0, samples], axis=1) - np.nanmedian(
            groups[1, samples], axis=1
        )
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    return {
        "depth": med[0] - med[1],
        "ci95": [lo, hi],
        "core_nights": int(counts[0]),
        "reference_nights": int(counts[1]),
        "reference_mag": med[1],
    }


def profile(params, p):
    level, a, b, ca, cb, wa, wb = params
    da = (p - ca + 0.5) % 1 - 0.5
    db = (p - cb) % 1 - 0.5
    return level + a * np.exp(-0.5 * (da / wa) ** 2) + b * np.exp(-0.5 * (db / wb) ** 2)


def fit_profile(p, mag, err):
    med = np.median(mag)
    lower = [med - 0.5, 0, 0, -0.12, -0.12, 0.008, 0.008]
    upper = [med + 0.5, 2, 2, 0.12, 0.12, 0.1, 0.1]
    fit = least_squares(
        lambda a: (profile(a, p) - mag) / np.hypot(err, 0.015),
        [med, 0.18, 0.18, 0, 0, 0.04, 0.04],
        bounds=(lower, upper),
        loss="soft_l1",
        f_scale=1.5,
    )
    return fit.x


def bounded_period_search(t, m, err, period):
    # Both half-period and orbital-period interpretations are covered via
    # a single box modeling eclipses repeating twice per orbit.
    from astropy.timeseries import BoxLeastSquares

    per = np.linspace(0.98, 1.02, 4001) * period / 2
    bls = BoxLeastSquares(t, -m, np.hypot(err, 0.02))
    res = bls.power(per, np.array([0.025, 0.05, 0.09]) * period, oversample=10)
    i = np.argmax(res.power)
    return {
        "eclipse_spacing_days": float(per[i]),
        "box_depth_mag": float(res.depth[i]),
        "duration_days": float(res.duration[i]),
    }


def run(raw, destination):
    catalog = {
        row.split()[0]: row.split()
        for row in (raw / "ecl.dat").read_text().splitlines()
    }
    results = {}
    rng = np.random.default_rng(111326)
    for number in STARS:
        star = "OGLE-SMC-ECL-" + number
        period, epoch = map(float, catalog[star][3:5])
        results[star] = {"period_days": period, "epoch_hjd": epoch, "bands": {}}
        for band in ["I", "V"]:
            times, mag, err = np.loadtxt(raw / "phot" / band / (star + ".dat")).T
            good = np.isfinite(times + mag + err) & (err > 0) & (err <= 0.1)
            p_all = phase(times, epoch, period)
            epochs = {}
            for i, era in enumerate(ERAS):
                use = good & (times >= EDGES[i]) & (times < EDGES[i + 1])
                t, p, m, e = [v[use] for v in (times, p_all, mag, err)]
                row = {"n": len(t), "nights": len(np.unique(np.floor(t)))}
                row["primary"] = depth_measure(t, p, m, 0, rng)
                row["secondary"] = depth_measure(t, p, m, 0.5, rng)
                if band == "I" and len(t) >= 80:
                    row["profile_parameters"] = fit_profile(p, m, e).tolist()
                    row["unphased_p95_minus_median_mag"] = float(
                        np.percentile(m, 95) - np.median(m)
                    )
                    # Sensitivity to rejecting observations with larger reported errors.
                    wider = (
                        (times >= EDGES[i])
                        & (times < EDGES[i + 1])
                        & (err > 0)
                        & (err <= 0.3)
                    )
                    row["wider_error_primary"] = depth_measure(
                        times[wider], p_all[wider], mag[wider], 0, rng
                    )
                    row["wider_error_secondary"] = depth_measure(
                        times[wider], p_all[wider], mag[wider], 0.5, rng
                    )
                epochs[era] = row
            results[star]["bands"][band] = epochs
            if band == "I" and number == "07730":
                last = good & (times >= EDGES[-2])
                early = epochs[ERAS[0]]["profile_parameters"]
                injected = profile(early, p_all[last]) - early[0]
                depth_search = {}
                for name, values in [
                    ("observed", mag[last]),
                    ("injected_early_profile", mag[last] + injected),
                ]:
                    depth_search[name] = bounded_period_search(
                        times[last], values, err[last], period
                    )
                results[star]["late_period_search"] = depth_search
    destination.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("output", type=Path)
    a = parser.parse_args()
    run(a.raw, a.output)
