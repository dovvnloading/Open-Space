"""Replay the HD 32564 discovery and its cross-instrument checks offline."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
EPOCH = 2456000.
PREDICTION_HASH = "441faa7b68b6635cda77cf9528a796b0cd5ccfb11afee7bef7579c541d1fff5a"


def design(t, periods):
    phase = 2 * np.pi * t[:, None] / np.asarray(periods)
    return np.stack((np.cos(phase), np.sin(phase)), axis=2).reshape(len(t), -1)


def remove(values, nuisance):
    coefficients = np.linalg.lstsq(nuisance, values, rcond=None)[0]
    return values - nuisance @ coefficients


def summaries(rows):
    bins = np.round(rows["BJD"])
    groups = [bins == day for day in np.unique(bins)]
    result = {}
    for key, values in rows.items():
        # RVBank uses -9999999 for missing measurements, including RHKp.
        valid = np.isfinite(values) & (values != -9999999)
        medians = []
        for group in groups:
            sample = values[group & valid]
            medians.append(np.median(sample) if sample.size else np.nan)
        result[key] = np.array(medians)
    return result


def fit_harps(t, y, seeds):
    nuisance = np.column_stack([np.ones(len(t)), (t + EPOCH > 2457165), t / 1000])

    def residual(periods):
        return remove(y, np.column_stack([design(t, periods), nuisance]))

    seeds = np.array(seeds)
    width = seeds**2 * .00015
    periods = least_squares(residual, seeds, bounds=(seeds-width, seeds+width),
                            xtol=1e-13, ftol=1e-13).x
    matrix = np.column_stack([design(t, periods), nuisance])
    coef = np.linalg.lstsq(matrix, y, rcond=None)[0]
    rms = np.std(y - matrix @ coef)
    return periods, coef, rms


def audit_signals(t, values, periods, nuisance):
    """Partial sinusoid fit at each fixed frequency; shared jitter is empirical."""
    matrix = np.column_stack([design(t, periods), nuisance])
    coef = np.linalg.lstsq(matrix, values, rcond=None)[0]
    residual = values - matrix @ coef
    output = []
    for i in range(len(periods)):
        restricted = np.delete(matrix, [2*i, 2*i+1], axis=1)
        r0 = remove(values, restricted)
        frac = 1 - np.sum(residual**2)/np.sum(r0**2)
        output.append({"K": round(float(np.hypot(*coef[2*i:2*i+2])), 6),
                       "partial_r2": round(float(frac), 4)})
    return output


def main():
    def numeric_columns(path, names):
        with path.open(newline="") as source:
            records = list(csv.DictReader(source))
        return {name: np.array([float(row[name]) for row in records]) for name in names}

    columns = ("BJD", "RV_mlc_nzp", "RV_drs_nzp", "CRX", "dLW", "Halpha",
               "RHKp", "FWHM_DRS", "BIS", "BERV")
    raw = numeric_columns(ROOT/"harps.csv", columns)
    h = summaries(raw)
    t, y = h["BJD"] - EPOCH, h["RV_mlc_nzp"]
    e = numeric_columns(ROOT/"espresso.csv", ("bjd", "rv", "rv_err", "bis",
                                             "contrast", "fwhm", "berv"))
    u, v = e["bjd"] - EPOCH, e["rv"] * 1000.
    base = np.column_stack([np.ones(len(u)), (u-u.mean())/1000])
    observed = remove(v, base)
    frozen = (ROOT/"harps_predictions.json").read_bytes()
    if hashlib.sha256(frozen).hexdigest() != PREDICTION_HASH:
        raise ValueError("Frozen HARPS predictions changed")
    rows = json.loads(frozen)
    aliases = []
    for row in rows:
        periods, coef, training_rms = fit_harps(t, y, row["seed"])
        np.testing.assert_allclose(periods, row["P"], rtol=0, atol=2e-6)
        np.testing.assert_allclose(coef, row["params"], rtol=0, atol=.001)
        # Only the four sinusoidal terms extrapolate. Instrument offset/trend
        # are nuisance terms fitted equally in the validation null and models.
        components = (design(u, row["P"]) * row["params"][:8]).reshape(-1, 4, 2).sum(2)
        predicted = remove(components.sum(1), base)
        miss = observed - predicted
        omission = [np.mean((miss + remove(wave, base))**2) - np.mean(miss**2)
                    for wave in components.T]
        aliases.append({"inner_days": round(row["P"][0], 6),
                        "harps_rms": round(float(training_rms), 5),
                        "espresso_r": round(float(np.corrcoef(observed,predicted)[0,1]), 5),
                        "espresso_rms": round(float(np.std(miss)), 5),
                        "omit_mse_increase": np.round(omission, 4).tolist()})
    chosen = min(rows, key=lambda row: row["rms"])
    periods = chosen["P"]
    ht = np.column_stack([np.ones(len(t)), t/1000, t+EPOCH>2457165])
    diagnostics = {}
    for name, values in h.items():
        if name == "BJD":
            continue
        valid = np.isfinite(values)
        diagnostics[name] = audit_signals(t[valid], values[valid], periods, ht[valid])
    # Independent coefficients for the four fixed-phase waves in ESPRESSO.
    unit_waves = (design(u, periods) * chosen["params"][:8]).reshape(-1, 4, 2).sum(2)
    inner_subset = {
        "shortest_only_rms": round(float(np.std(remove(v - unit_waves[:,0], base))), 5),
        "two_shortest_rms": round(float(np.std(remove(v - unit_waves[:,:2].sum(axis=1), base))), 5),
    }
    unit_waves /= np.array(chosen["K"])
    matrix = np.column_stack([unit_waves, base])
    coef = np.linalg.lstsq(matrix, v, rcond=None)[0]
    residual = v - matrix @ coef
    covariance = np.linalg.inv(matrix.T @ matrix) * (residual @ residual)/(len(v)-matrix.shape[1])
    espresso_amplitude = {"K_m_s": np.round(coef[:4], 4).tolist(),
                          "ols_se_m_s": np.round(np.sqrt(covariance.diagonal()[:4]), 4).tolist()}
    espresso_indicators = {}
    for name in ("bis", "contrast", "fwhm", "berv"):
        response = remove(e[name], base)
        projection = remove(unit_waves @ np.array(chosen["K"]), base)
        espresso_indicators[name] = round(float(np.corrcoef(response, projection)[0,1]), 4)
    rng = np.random.default_rng(32564)
    pred = design(u, periods) @ chosen["params"][:8]
    # Remove group means, then shuffle predictions within observing seasons.
    # This preserves actual seasonal means, sampling and measured RVs.
    years = np.floor((u + EPOCH - 2460175)/365.25).astype(int)
    groups = [np.flatnonzero(years == year) for year in np.unique(years)]
    group_means = np.column_stack([years == year for year in np.unique(years)])
    measured = remove(v, group_means)
    projected = remove(pred, group_means)
    score = measured @ projected
    n_draws = 100000
    count = 0
    for _ in range(n_draws//1000):
        products = np.zeros(1000)
        for idx in groups:
            shuffled = np.argsort(rng.random((1000, len(idx))), axis=1)
            products += (measured[idx] * projected[idx][shuffled]).sum(axis=1)
        count += int(np.sum(products >= score))
    return {"harps_exposures":len(raw["BJD"]),"harps_rounded_bjd_bins":len(t),"espresso_exposures":len(u),
            "espresso_null_rms":round(float(np.std(observed)),5),
            "aliases":aliases,"harps_control_fits":diagnostics,
            "espresso_frozen_inner_subsets":inner_subset,
            "espresso_fixed_phase":espresso_amplitude,"espresso_indicator_r":espresso_indicators,
            "espresso_control":{"season_sizes":list(map(len,groups)),
                                "correlation":round(float(np.corrcoef(measured,projected)[0,1]),5),
                                "shuffles":n_draws,"exceedances":count}}


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, allow_nan=False))
