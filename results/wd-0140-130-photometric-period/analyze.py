"""Reproduce the period measurement using Python, NumPy, SciPy and Astropy."""

import csv
import gzip
import json
from pathlib import Path

import numpy as np
from astropy.timeseries import LombScargle
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
REF = 2458000.0
SIDEREAL = 1.00273790935


def read(name):
    with gzip.open(HERE / name, "rt") as stream:
        return list(csv.DictReader(stream))


def clean_ztf(rows):
    groups = {}
    for row in rows:
        mag, err = float(row["mag"]), float(row["magerr"])
        if row["catflags"] == "0" and err < 0.3 and mag < float(row["limitmag"]) - 0.2:
            groups.setdefault(row["oid"], []).append(row)
    retained = {}
    for detections in groups.values():
        if len(detections) < 12:
            continue
        mags = np.array([float(row["mag"]) for row in detections])
        median = np.median(mags)
        scale = max(0.025, 1.4826 * np.median(abs(mags - median)))
        for row, mag in zip(detections, mags):
            if abs(mag - median) > 5 * scale:
                continue
            key = row["filtercode"], row["expid"]
            measurement = (
                float(row["hjd"]) - REF,
                mag - median,
                np.hypot(float(row["magerr"]), 0.015),
                row["filtercode"],
            )
            if key not in retained or measurement[2] < retained[key][2]:
                retained[key] = measurement
    return tuple(map(np.array, zip(*retained.values())))


def design(t, band, freq):
    angle = 2 * np.pi * t * freq
    return np.column_stack(
        [
            v * (band == k)
            for k in np.unique(band)
            for v in (np.ones(len(t)), np.cos(angle), np.sin(angle))
        ]
    )


def fit(t, y, e, b, freq):
    X = design(t, b, freq)
    beta = np.linalg.lstsq(X / e[:, None], y / e, rcond=None)[0]
    pred = X @ beta
    baseline = np.zeros(len(y))
    for k in np.unique(b):
        mask = b == k
        baseline[mask] = np.average(y[mask], weights=1 / e[mask] ** 2)
    gain = ((y - baseline) ** 2 - (y - pred) ** 2) / e**2
    return float(sum(gain)), beta, pred


def subset(data, mask):
    return tuple(col[mask] for col in data)


def fine(data, freq):
    step = 0.5 / np.ptp(data[0])
    return minimize_scalar(
        lambda f: -fit(*data, f)[0],
        bounds=(freq - step, freq + step),
        method="bounded",
        options={"xatol": 1e-12},
    ).x


def bootstrap_alias(t, y, e, b, primary, alternate):
    # Resample observed nights of fixed fitted models. This quantifies model preference,
    # conditional on these two aliases and their fitted coefficients.
    p = fit(t, y, e, b, primary)[2]
    alt = fit(t, y, e, b, alternate)[2]
    _, inv = np.unique(np.floor(t + 0.5).astype(int), return_inverse=True)
    per_night = np.bincount(inv, weights=((y - alt) ** 2 - (y - p) ** 2) / e**2)
    rng = np.random.default_rng(764132)
    draws = np.concatenate(
        [rng.choice(per_night, (100, len(per_night))).sum(axis=1) for _ in range(50)]
    )
    return {
        "nights": len(per_night),
        "delta_chi2": sum(per_night),
        "interval95": np.quantile(draws, [0.025, 0.975]).tolist(),
        "fraction_prefer_primary": float(np.mean(draws > 0)),
    }


def run():
    raw = read("ztf.csv.gz")
    data = clean_ztf(raw)
    gr = subset(data, data[3] != "zi")
    t, y, e, b = gr
    split = np.quantile(t, 0.65)
    early, late = subset(gr, t <= split), subset(gr, t > split)
    # Equal frequency grid with five samples per independent Fourier spacing.
    fe = np.arange(0.15, 65, 1 / (5 * np.ptp(early[0])))
    periodogram = np.zeros(len(fe))
    for k in ["zg", "zr"]:
        tk, yk, ek, _ = subset(early, early[3] == k)
        periodogram += LombScargle(tk, yk, ek, normalization="psd").power(
            fe, method="fast", assume_regular_frequency=True
        )
    blind_peak = fe[periodogram.argmax()]
    # Compare the three neighboring sidereal aliases of the early maximum.
    family = [
        blind_peak + j * SIDEREAL
        for j in (-2, -1, 0, 1)
        if blind_peak + j * SIDEREAL > 0.15
    ]
    combined = [fine(gr, candidate) for candidate in family]
    f = max(combined, key=lambda candidate: fit(*gr, candidate)[0])
    f_early = fine(early, f)
    gain_train, beta, _ = fit(*early, f_early)
    predicted = design(late[0], late[3], f_early) @ beta
    null = np.array([beta[0] if k == "zg" else beta[3] for k in late[3]])
    statistic = ((late[1] - null) ** 2 - (late[1] - predicted) ** 2) / late[2] ** 2
    _, inv = np.unique(np.floor(late[0] + 0.5).astype(int), return_inverse=True)
    # Null symmetry test. Preserve all observations of a night together. Remove the
    # model-energy term first so the permutation flips residual flux, not a constant.
    linear = np.bincount(
        inv, weights=2 * (late[1] - null) * (predicted - null) / late[2] ** 2
    )
    rng = np.random.default_rng(808122)
    exceed = 0
    for _ in range(200):
        sums = rng.choice([-1, 1], (500, len(linear))) @ linear
        exceed += np.sum(sums >= linear.sum())
    aliases = {}
    for label, candidate in [
        ("signal", f),
        ("plus_sidereal", f + SIDEREAL),
        ("minus_sidereal", f - SIDEREAL),
        ("plus_solar", f + 1),
        ("minus_solar", f - 1),
        ("half", f / 2),
        ("double", f * 2),
    ]:
        ff = fine(gr, candidate)
        aliases[label] = {
            "frequency": ff,
            "hours": 24 / ff,
            "gain": fit(*gr, ff)[0],
            "by_band": {
                k: fit(*subset(data, data[3] == k), ff)[0] for k in ["zg", "zr", "zi"]
            },
        }
    # Local curvature only; separately fit each band to expose window or calibration differences.
    h = 1e-6
    center = fit(*gr, f)[0]
    sigma = np.sqrt(-2 * h * h / (fit(*gr, f - h)[0] - 2 * center + fit(*gr, f + h)[0]))
    by_band = {}
    for k in ["zg", "zr", "zi"]:
        sample = subset(data, data[3] == k)
        _, beta_k, _ = fit(*sample, f)
        by_band[k] = {
            "n": len(sample[0]),
            "fixed_semimagnitude": np.hypot(*beta_k[1:]),
            "fixed_phase_rad": np.arctan2(beta_k[2], beta_k[1]),
            "local_frequency": fine(sample, f),
        }
    halves = {}
    for name, s in [("early", early), ("late", late)]:
        # Free local peak for signal / alias in each interval; no fixed-phase constraint.
        ref = fit(*s, fine(s, f))[0]
        sid = fit(*s, fine(s, f + SIDEREAL))[0]
        halves[name] = {
            "n": len(s[0]),
            "local_signal_gain": ref,
            "local_plus_sidereal_gain": sid,
        }
    ps = []
    for row in read("panstarrs.csv.gz"):
        flux, err = float(row["psfFlux"]), float(row["psfFluxErr"])
        if (
            float(row["psfQfPerfect"]) > 0.85
            and flux / err > 10
            and not int(row["infoFlag3"]) & (0x01000000 | 0x00002000 | 0x00000004)
        ):
            ps.append(
                (
                    float(row["obsTime"]),
                    -2.5 * np.log10(flux),
                    np.hypot(1.085736 * err / flux, 0.01),
                    row["filterID"],
                )
            )
    ps = tuple(map(np.array, zip(*ps)))
    return {
        "raw_rows": len(raw),
        "total_gr": len(t),
        "baseline_days": np.ptp(t),
        "split_hjd": split + REF,
        "blind_early_highest_frequency": float(blind_peak),
        "early_alias_family": family,
        "train_frequency": f_early,
        "train_gain": gain_train,
        "holdout_fixed_gain": sum(statistic),
        "holdout_refit_gain": fit(*late, f_early)[0],
        "holdout_nights": len(linear),
        "holdout_randomization": {
            "trials": 100000,
            "as_or_more_extreme": int(exceed),
            "plus_one_p": (1 + int(exceed)) / 100001,
        },
        "frequency": f,
        "hours": 24 / f,
        "formal_hours_error": 24 * sigma / f**2,
        "aliases": aliases,
        "night_bootstrap_plus_sidereal": bootstrap_alias(
            *gr, f, aliases["plus_sidereal"]["frequency"]
        ),
        "by_band": by_band,
        "halves": halves,
        "ps1_n": len(ps[0]),
        "ps1_band_free_gain": {
            k: fit(*ps, v["frequency"])[0] for k, v in aliases.items()
        },
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
