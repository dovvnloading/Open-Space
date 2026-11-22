"""Targeted published-orbit comparisons on the survey's original RVBank nights.

The blind period was selected on these same data. BIC differences and the ordered
split are diagnostics, with no survey-wide significance assigned to the added orbit.
"""
import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from runtime import require_survey_workflow
from search import nightly, orbit_basis

TAU = 2 * np.pi


@dataclass(frozen=True)
class Planet:
    period: float
    sigma: float
    eccentricity: float


# The first planet is the published counterpart under test. HD 216770 c covers
# a small fraction of its reported orbit here; a quadratic is our default, and
# the second baseline uses the thesis's nominal outer orbit instead.
SYSTEMS = {
    "GJ550": (1, [Planet(960.41, 6.33, .035)], "https://www.aanda.org/articles/aa/pdf/2018/07/aa32791-18.pdf"),
    "HD216770": (1, [Planet(118.45, .55, .37)],
                 "https://doi.org/10.1051/0004-6361:20034250 | https://archive-ouverte.unige.ch/unige:141226"),
    "GJ634.1": (1, [Planet(1096.27, 27.24, 0)], "https://arxiv.org/abs/1109.2497"),
    "GJ496.1": (2, [Planet(663.2, 8., .14), Planet(1818., 25., .20)],
                "https://doi.org/10.1051/0004-6361/201424965"),
}


@dataclass(frozen=True)
class Orbit:
    seed: Planet
    low: float
    high: float
    circular: bool = False


def orbit_family(planets, blind, mode):
    orbits = []
    for idx, planet in enumerate(planets):
        p = planet.period
        low = high = p
        if mode == "blind replacement" and idx == 0:
            p = low = high = blind
        elif mode == "alternative published period" and idx == 0:
            p = low = high = .318 * 365.25
        elif mode.startswith("published range"):
            low, high = p - 3 * planet.sigma, p + 3 * planet.sigma
        elif mode.startswith("broad range"):
            low, high = p * .85, p * 1.15
        orbits.append(Orbit(Planet(p, planet.sigma, planet.eccentricity), low, high))
    if "+ extra" in mode:
        # Fixed blind frequency, with both sine and cosine fitted. This residual
        # test uses two parameters and targets a circular component specifically.
        orbits.append(Orbit(Planet(blind, 0, 0), blind, blind, circular=True))
    return orbits


def build_model(epochs, mask, family, baseline, rng):
    """Optimize orbital geometry and solve offsets, drift and amplitudes linearly."""
    times = epochs.bjd.to_numpy()
    eras = epochs.era.to_numpy()
    y, err = epochs.rv.to_numpy(), epochs.err.to_numpy()
    t = times - np.mean(times[mask])
    span = np.ptp(t[mask])
    base = [np.ones(len(t)), t / span]
    if len(np.unique(eras[mask])) > 1:
        base.append((eras == "POST").astype(float))
    if baseline == "quadratic outer arc":
        base.append((t / span)**2)
    elif baseline == "nominal outer Keplerian":
        # Rickman 2020 Table 4.10: 102.56 yr, e=.17, T0=BJD 2476121.
        # The two fitted coefficients absorb K and argument of periastron.
        base.extend(orbit_basis(times - 2476121., 1 / (102.56 * 365.25), .17, 0).T)

    lower, upper, initial = [], [], []
    for orb in family:
        if orb.low < orb.high:
            lower.append(np.log(orb.low / orb.seed.period))
            upper.append(np.log(orb.high / orb.seed.period))
            initial.append(0.)
        if not orb.circular:
            lower.extend([0., -2*TAU])
            upper.extend([.80, 2*TAU])
            initial.extend([max(.03, orb.seed.eccentricity), 0.])

    def evaluate(par, include_prediction=False):
        j = 0
        columns = list(base)
        periods, eccentricities = [], []
        for orb in family:
            p = orb.seed.period
            if orb.low < orb.high:
                p *= np.exp(par[j])
                j += 1
            e, phase = 0., 0.
            if not orb.circular:
                e, phase = par[j:j+2]
                j += 2
            columns.extend(orbit_basis(t, 1 / p, e, phase).T)
            periods.append(p)
            eccentricities.append(e)
        design = np.column_stack(columns)
        a = design[mask] / err[mask, None]
        coef = np.linalg.lstsq(a, y[mask] / err[mask], rcond=None)[0]
        prediction = design @ coef
        res = (prediction[mask] - y[mask]) / err[mask]
        if include_prediction:
            return prediction, periods, eccentricities, coef, a, res
        return res

    best = None
    lower, upper = np.asarray(lower), np.asarray(upper)
    count = 12 if len(family) >= 2 else 7
    for start in range(count if initial else 1):
        par = np.asarray(initial, dtype=float)
        j = 0
        for orb in family:
            if orb.low < orb.high:
                if start >= 2:
                    par[j] = rng.uniform(lower[j], upper[j])
                j += 1
            if not orb.circular:
                if start > 0:
                    par[j:j+2] = [rng.uniform(.04, .60), rng.uniform(-np.pi, np.pi)]
                j += 2
        if len(par):
            fit = least_squares(evaluate, par, bounds=(lower, upper), max_nfev=160,
                                xtol=3e-8, ftol=3e-8, gtol=3e-8)
            par = fit.x
            converged = fit.success
        else:
            converged = True
        residual = evaluate(par)
        rss = float(residual @ residual)
        if best is None or rss < best[0]:
            best = rss, par, converged
    rss, params, converged = best
    pred, periods, eccs, coeff, a, _ = evaluate(params, True)
    n = mask.sum()
    npars = a.shape[1] + len(params) + 1  # one common noise scale
    amps = np.hypot(coeff[len(base)::2], coeff[len(base)+1::2])
    return dict(
        prediction=pred, rss=rss, bic=n*np.log(max(rss/n, 1e-30)) + npars*np.log(n),
        npars=npars, weighted_rms_ms=np.sqrt(rss / np.sum(1 / err[mask]**2)),
        period=periods, eccentricity=eccs, amplitude=amps,
        linear_condition_number=np.linalg.cond(a / np.linalg.norm(a, axis=0)),
        optimizer_converged=converged,
    )


def score_tail(epochs, fit, training):
    # The blind frequencies were selected on the entire release, so this
    # chronological check cannot provide an independent false-alarm probability.
    test = ~training & epochs.era.isin(epochs.era[training]).to_numpy()
    if test.sum() < 12:
        return int(test.sum()), np.nan
    y = epochs.rv.to_numpy() - fit["prediction"]
    err = epochs.err.to_numpy()
    for era in set(epochs.era[test]):
        rows = test & epochs.era.eq(era).to_numpy()
        y[rows] -= np.average(y[rows], weights=1 / err[rows]**2)
    return int(test.sum()), np.sqrt(np.average(y[test]**2, weights=1 / err[test]**2))


def values(nums):
    return " | ".join(f"{x:.8g}" for x in nums)


def main(args):
    require_survey_workflow()
    args.out.mkdir(parents=True, exist_ok=True)
    survivors = pd.read_csv(args.evidence / "rv-survivors.csv").set_index(["star","signal_index"])
    holdout = pd.read_csv(args.evidence / "holdout.csv").set_index(["star","signal_index"])
    target_chunks = []
    cols = ["star","bjd","rv","err","era"]
    for chunk in pd.read_csv(args.primary, usecols=cols, chunksize=100000):
        target_chunks.append(chunk[chunk.star.isin(SYSTEMS)])
    groups = pd.concat(target_chunks).groupby("star")
    rows = []
    for star, (idx, planets, source) in SYSTEMS.items():
        blind = survivors.loc[(star, idx), "period_d"]
        later = holdout.loc[(star, idx)]
        epochs = nightly(groups.get_group(star))
        size = len(epochs)
        training = np.arange(size) < int(.75 * size)
        all_rows = np.ones(size, dtype=bool)
        baseline_modes = ["linear drift"] if star != "HD216770" else [
            "quadratic outer arc", "nominal outer Keplerian"]
        span = np.ptp(epochs.bjd)
        for baseline in baseline_modes:
            modes = ["nuisance only", "published fixed", "blind replacement", "published range",
                     "published range + extra", "broad range", "broad range + extra"]
            if star == "HD216770":
                modes.insert(3, "alternative published period")
            model_rows = {}
            for mode in modes:
                orbits = [] if mode=="nuisance only" else orbit_family(planets, blind, mode)
                fits = []
                for mask in [all_rows, training]:
                    fits.append(build_model(epochs, mask, orbits, baseline, np.random.default_rng(46811)))
                full, train = fits
                n_test, wrms = score_tail(epochs, train, training)
                row = dict(star=star, signal_index=idx, blind_period_d=blind, baseline=baseline,
                           model=mode, nights=size, baseline_d=span,
                           blind_to_published_rayleigh=abs(1/blind-1/planets[0].period)*span,
                           tested_periods_d=values(full["period"]), eccentricities=values(full["eccentricity"]),
                           semiamplitudes_ms=values(full["amplitude"]), n_model_parameters=full["npars"],
                           bic=full["bic"], fitted_weighted_rms_ms=full["weighted_rms_ms"],
                           optimizer_converged=full["optimizer_converged"],
                           linear_condition_number=full["linear_condition_number"],
                           train_nights=int(training.sum()), scored_tail_nights=n_test,
                           tail_shape_wrms_ms=wrms, tail_optimizer_converged=train["optimizer_converged"],
                           independent_later_nights=int(later.n_later_nights), independent_outcome=later.outcome,
                           source=source)
                if "+ extra" in mode:
                    row["extra_to_fitted_nearest_rayleigh"] = min(
                        abs(1/blind - 1/p)*span for p in full["period"][:-1])
                model_rows[mode] = row
            for mode, row in model_rows.items():
                row["delta_bic_vs_published_fixed_positive_favors_row"] = model_rows["published fixed"]["bic"] - row["bic"]
                if "+ extra" in mode:
                    base = model_rows[mode.split(" +")[0]]
                    row["delta_bic_added_positive_favors_two"] = base["bic"] - row["bic"]
                    row["tail_wrms_improvement_added_ms"] = base["tail_shape_wrms_ms"] - row["tail_shape_wrms_ms"]
                rows.append(row)
        print(f"{star}: {size} nights and {len(baseline_modes)} nuisance baseline(s)", flush=True)
    pd.DataFrame(rows).to_csv(args.out / "prior-orbit-comparison.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
