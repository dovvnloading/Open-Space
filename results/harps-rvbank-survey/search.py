"""Dispatch-only HARPS time-series search. All velocities are in m/s and times in days."""
import argparse
import csv
import hashlib
import json
import os
from contextlib import ExitStack
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from scipy.stats import f as f_dist, spearmanr

INDICATORS = ("fwhm", "bis", "crx", "dlw", "halpha", "sindex", "rhkp")
REQUIRED = INDICATORS[:6]
MIN_NIGHTS = 12
MIN_SPAN = 30.0
MIN_PERIOD = 1.5
RESAMPLES = 255
MAX_SIGNALS = 3
TAU = 2 * np.pi
STAR_COLUMNS = ("star n_epochs n_nights baseline_d n_pre n_post eligibility "
                "slope_ms_yr slope_err_ms_yr trend_bjd n_freq n_scans n_blocks "
                "period_min_d period_max_d error").split()
SIGNAL_COLUMNS = ("star signal_index period_d k_ms ecc omega_rad tperi_bjd grid_period_d "
                  "fap fap_gaussian_grid gaussian_trials_screen peak_power "
                  "delta_bic phase_gap offset_post_minus_pre_ms bootstrap_exceedances "
                  "n_bootstrap n_blocks window_peaks_d alias_periods_d alias_fits_json activity_json "
                  "amplitude_checks_json rv_status flags").split()
INJECTION_COLUMNS = ("star period_d k_ms eccentricity trials recovered "
                     "gls_recovery_fraction threshold_power").split()


def nightly(df):
    """Use the integer Julian day, which changes during La Silla daylight."""
    g = df.copy()
    for col in ["bjd", "rv", "err", *INDICATORS]:
        g[col] = pd.to_numeric(g.get(col, np.nan), errors="coerce")
    valid = np.isfinite(g[["bjd", "rv", "err"]]).all(axis=1) & (g.err > 0)
    g = g.loc[valid].copy()
    if g.empty:
        return g
    if "era" not in g or not g.era.isin(["PRE", "POST"]).all():
        raise ValueError("every valid epoch needs an explicit PRE or POST fibre era")
    g["day"] = np.floor(g.bjd).astype(int)
    w = 1 / g.err**2
    keys = [g.day, g.era]
    wg = w.groupby(keys).sum()
    data = {col: (g[col] * w).groupby(keys).sum() / wg for col in ("bjd", "rv")}
    # Add a 1 m/s nightly floor; use the larger of formal and internal scatter SEM.
    count = g.rv.groupby(keys).count()
    sem = g.rv.groupby(keys).std().fillna(0) / np.sqrt(count)
    data["err"] = np.sqrt(1 + np.maximum(sem**2, 1 / wg))
    for key in INDICATORS:
        data[key] = g[key].replace([np.inf, -np.inf], np.nan).groupby(keys).median()
    return pd.DataFrame(data).reset_index().sort_values("bjd").reset_index(drop=True)


def nuisance(t, era):
    cols = [np.ones(len(t)), (t - np.mean(t)) / 365.25]
    if np.unique(era).size > 1:
        cols.append((era == "POST").astype(float))
    return np.column_stack(cols)


def blocks_for(t, era):
    block = np.floor((t - t.min()) / 14).astype(int)
    return np.unique(np.stack([block, era == "POST"], axis=1), axis=0,
                     return_inverse=True)[1].ravel()


def trend(t, y, err, design, blocks):
    a = design / err[:, None]
    coef = np.linalg.lstsq(a, y / err, rcond=None)[0]
    res = (y - design @ coef) / err
    inv = np.linalg.pinv(a.T @ a)
    leverage = np.sum((a @ inv) * a, axis=1)
    scores = a * (res / np.maximum(1 - leverage, .01))[:, None]
    score_by_block = np.zeros((blocks.max() + 1, len(coef)))
    np.add.at(score_by_block, blocks, scores)
    cov = (score_by_block @ inv).T @ (score_by_block @ inv)
    return coef, np.sqrt(np.maximum(np.diag(inv), np.diag(cov)))


class Grid:
    """Weighted sine/cosine regression after projecting out every null-model term."""

    def __init__(self, t, err, design, frequencies):
        self.t = t - t.min()
        self.err = err
        self.f = frequencies
        self.q = np.linalg.qr(design / err[:, None], mode="reduced")[0]

    def project(self, series):
        z = np.atleast_2d(series) / self.err
        return z - (z @ self.q) @ self.q.T

    def scan(self, series, return_first=False):
        z = self.project(series)
        raw = np.atleast_2d(series) / self.err
        # Exact offset/drift-only data leave roundoff after QR projection.
        # Treat that as zero before normalizing by the residual variance.
        roundoff_limit = np.linalg.norm(raw, axis=1) * (128 * np.finfo(float).eps)
        z[np.linalg.norm(z, axis=1) <= roundoff_limit] = 0
        norm = np.maximum(np.einsum("ij,ij->i", z, z), np.finfo(float).tiny)
        maximum = np.full(len(z), -1.)
        best = np.zeros(len(z), dtype=int)
        first = []
        # Bound temporary memory even on the densest time series.
        width = max(32, min(2048, 2_000_000 // max(len(self.t), len(z))))
        for start in range(0, len(self.f), width):
            ff = self.f[start:start + width]
            ang = TAU * self.t[:, None] * ff
            s, c = np.sin(ang) / self.err[:, None], np.cos(ang) / self.err[:, None]
            s -= self.q @ (self.q.T @ s)
            c -= self.q @ (self.q.T @ c)
            ss = np.einsum("ij,ij->j", s, s)
            cc = np.einsum("ij,ij->j", c, c)
            sc = np.einsum("ij,ij->j", s, c)
            det = ss * cc - sc**2
            ok = (det > 1e-10 * (ss + cc)**2) & (det > 0)
            ys, yc = z @ s, z @ c
            gain = (cc * ys**2 + ss * yc**2 - 2 * sc * ys * yc)
            power = np.clip(gain / np.where(ok, det, 1) / norm[:, None], 0, 1)
            power[:, ~ok] = 0
            idx = power.argmax(axis=1)
            peak = power[np.arange(len(power)), idx]
            improve = peak > maximum
            best[improve] = start + idx[improve]
            maximum = np.maximum(maximum, peak)
            if return_first:
                first.extend(power[0])
        return maximum, best, np.asarray(first)


def period_grid(span):
    step = 1 / (5 * span)
    return np.arange(1 / span, 1 / MIN_PERIOD + step * 1e-8, step)


def window_peaks(t, err, span):
    freq = np.arange(2 / span, 2.05, 1 / (5 * span))
    w = (1 / err**2) / np.sum(1 / err**2)
    power = []
    for i in range(0, len(freq), max(32, min(2048, 1_000_000 // len(t)))):
        sub = freq[i:i + max(32, min(2048, 1_000_000 // len(t)))]
        ang = TAU * sub[:, None] * (t - t.min())
        power.extend((np.cos(ang) @ w)**2 + (np.sin(ang) @ w)**2)
    power = np.asarray(power)
    indices = find_peaks(power, height=.12, distance=4)[0]
    return freq[indices[np.argsort(power[indices])[-16:]]]


def orbit_basis(t, f, ecc, m0):
    m = (TAU * f * t + m0 + np.pi) % TAU - np.pi
    e_anom = m + .85 * ecc * np.sign(np.sin(m))
    for _ in range(12):
        change = (e_anom - ecc * np.sin(e_anom) - m) / (1 - ecc * np.cos(e_anom))
        e_anom -= change
        if np.max(np.abs(change)) < 1e-13:
            break
    den = 1 - ecc * np.cos(e_anom)
    return np.column_stack([(np.cos(e_anom) - ecc) / den + ecc,
                            np.sqrt(1 - ecc**2) * np.sin(e_anom) / den])


def fit_orbit(t, y, err, design, freq, span):
    times = t - t.min()
    # Search in Rayleigh-frequency units to keep nonlinear parameters comparable.
    lo = max(-.8, (1 / span - freq) * span)
    hi = min(.8, (1 / MIN_PERIOD - freq) * span)

    def evaluate(pars):
        df, ecc, m0 = pars
        bases = orbit_basis(times, freq + df / span, ecc, m0)
        a = np.column_stack([design, bases])
        coef = np.linalg.lstsq(a / err[:, None], y / err, rcond=None)[0]
        return coef, (y - a @ coef) / err, bases

    solutions = []
    for ecc, phase in [(0.03, 0.), (.35, 2.), (.55, 4.)]:
        sol = least_squares(lambda p: evaluate(p)[1],
                            [np.clip(0., lo + 1e-8, hi - 1e-8), ecc, phase],
                            bounds=([lo, 0, -TAU * 2], [hi, .85, TAU * 2]),
                            max_nfev=90)
        solutions.append(sol)
    sol = min(solutions, key=lambda fit: fit.cost)
    df, ecc, m0 = sol.x
    coef, res, bases = evaluate(sol.x)
    period = 1 / (freq + df / span)
    k = np.hypot(*coef[-2:])
    omega = np.arctan2(-coef[-1], coef[-2])
    tperi = t.min() - (m0 % TAU) * period / TAU
    return period, k, ecc, omega, tperi, coef, res, bases, sol.success


def activity(night, y, freq, span):
    out, flags = {}, []
    for name in INDICATORS:
        good = np.isfinite(night[name].to_numpy())
        count = int(good.sum())
        data = {"n": count}
        out[name] = data
        values = night[name].to_numpy()[good]
        if count < 10 or np.ptp(values) == 0:
            data["status"] = "missing_or_insufficient"
            if name in REQUIRED:
                flags.append("missing_" + name)
            continue
        t = night.bjd.to_numpy()[good]
        era = night.era.to_numpy()[good]
        design = nuisance(t, era)
        a = Grid(t, np.ones(count), design, np.array([freq, 2 * freq, freq / 2]))
        val_res = a.project(values)[0]
        rv_res = a.project(y[good])[0]
        rho, rho_p = spearmanr(val_res, rv_res)
        _, _, power = a.scan(values, return_first=True)
        dof = count - design.shape[1] - 2
        p = f_dist.sf(power * dof / (2 * np.maximum(1 - power, 1e-14)), 2, dof)
        data.update(rho=float(rho), rho_p=float(rho_p), harmonic_p=[float(v) for v in p])
        data["status"] = "tested"
        if (abs(rho) >= .4 and rho_p < .01 / len(INDICATORS)) or min(p) < .01 / (3 * len(INDICATORS)):
            flags.append("activity_" + name)
        # Search for a rotation proxy without assigning every activity peak to rotation.
        span_act = np.ptp(t)
        frot = np.arange(max(1 / 200., 1 / span_act), 1 / MIN_PERIOD, 1 / (5 * span_act))
        if len(frot):
            a = Grid(t, np.ones(count), design, frot)
            high, idx, _ = a.scan(values)
            actf = frot[idx[0]]
            peak = high[0]
            pgrid = min(1., len(frot) * f_dist.sf(peak * dof / (2 * max(1 - peak, 1e-14)), 2, dof))
            data.update(proxy_period_d=1 / actf, proxy_grid_fap=pgrid)
            if pgrid < .01 / len(INDICATORS) and min(abs(freq - actf * h) for h in (.5, 1., 2., 3.)) <= 1 / span:
                flags.append("rotation_proxy_" + name)
    return out, flags


def aliases(grid, powers, windows, span, fitted_freq):
    high = powers.max()
    idx = find_peaks(powers)[0]
    idx = np.unique(np.r_[idx, 0, len(powers) - 1])
    idx = idx[(powers[idx] >= .8 * high) & (abs(grid[idx] - fitted_freq) > 1 / span)]
    alternatives = grid[idx[np.argsort(powers[idx])[::-1]]]
    competing = []
    for alt in alternatives:
        if np.any(np.abs(windows - abs(alt - fitted_freq)) < 1 / span) or np.any(abs(windows - (alt + fitted_freq)) < 1 / span):
            competing.append(float(1 / alt))
    return competing, bool(np.any(abs(windows - fitted_freq) < 1 / span))


def amplitude_checks(t, y, err, era, shape):
    """Independent fixed-shape amplitudes in chronological halves and, if sampled, eras."""
    checks = {}
    flags = []
    partitions = {"time": t > np.median(t)}
    if len(set(era)) == 2:
        partitions["era"] = era == "POST"
    for name, mask in partitions.items():
        amps = []
        for g in [~mask, mask]:
            if sum(g) < 5:
                continue
            a = np.column_stack([nuisance(t[g], era[g]), shape[g]])
            aw = a / err[g, None]
            if np.linalg.matrix_rank(aw) != a.shape[1] or np.linalg.cond(aw) > 1e7:
                continue
            coef, std = trend(t[g], y[g], err[g], a, blocks_for(t[g], era[g]))
            amps.append([float(coef[-1]), float(std[-1])])
        checks[name] = amps
        if len(amps) != 2:
            flags.append(name + "_split_unresolved")
        elif abs(amps[0][0] - amps[1][0]) > 3 * np.hypot(amps[0][1], amps[1][1]):
            flags.append(name + "_split_inconsistent")
    return checks, flags


def inject(star, grid, residual, threshold, rng):
    span = np.ptp(grid.t)
    series, cells = [], []
    for period in [3., 10., 30., 100., 365.25, 1000., 3000.]:
        if period > span:
            continue
        for ecc in [0., .35]:
            # Stratify mean anomaly with a random phase within each eighth of a cycle.
            phases = TAU * (np.arange(8) + rng.random(8)) / 8
            waves = []
            for phase in phases:
                omega = rng.uniform(0, TAU)
                waves.append(orbit_basis(grid.t, 1 / period, ecc, phase) @ [np.cos(omega), -np.sin(omega)])
            waves = np.asarray(waves)
            for k in [1., 3., 10., 30., 100., 300.]:
                series.extend(residual + k * waves)
                cells.append((period, k, ecc))
    powers, idx, _ = grid.scan(series)
    rows = []
    for i, (period, k, ecc) in enumerate(cells):
        ix = slice(8 * i, 8 * (i + 1))
        tol = min(.1 / period, 1 / span)
        count = int(np.count_nonzero((powers[ix] > threshold) & (abs(grid.f[idx[ix]] - 1 / period) <= tol)))
        rows.append(dict(star=star, period_d=period, k_ms=k, eccentricity=ecc,
                         trials=8, recovered=count, gls_recovery_fraction=count / 8,
                         threshold_power=threshold))
    return rows


def search(star, df, star_count):
    rng = np.random.default_rng(int.from_bytes(hashlib.sha256(star.encode()).digest()[:8], "little"))
    n = nightly(df)
    t = n.bjd.to_numpy()
    span = float(np.ptp(t)) if len(t) else 0.
    summary = dict(star=star, n_epochs=len(df), n_nights=len(t), baseline_d=span,
                   n_pre=int((n.era == "PRE").sum()), n_post=int((n.era == "POST").sum()))
    if len(t) < MIN_NIGHTS or span < MIN_SPAN:
        summary["eligibility"] = "insufficient_nights" if len(t) < MIN_NIGHTS else "insufficient_span"
        return summary, [], []
    era, y, err = n.era.to_numpy(), n.rv.to_numpy(), n.err.to_numpy()
    design = nuisance(t, era)
    blocks = blocks_for(t, era)
    n_blocks = int(blocks.max() + 1)
    base = design.copy()
    freq = period_grid(span)
    summary.update(eligibility="searched", n_blocks=n_blocks, period_min_d=1 / freq[-1],
                   period_max_d=1 / freq[0], n_freq=len(freq), trend_bjd=np.mean(t))
    windows = window_peaks(t, err, span)
    signals = []
    for seq in range(1, MAX_SIGNALS + 1):
        grid = Grid(t, err, design, freq)
        residual = grid.project(y)[0] * err
        signs = rng.integers(0, 2, (RESAMPLES, n_blocks)) * 2 - 1
        high, idx, powers = grid.scan(np.vstack([y, residual * signs[:, blocks]]), True)
        peak = high[0]
        exc = int((high[1:] >= peak - 1e-12).sum())
        fap = (1 + exc) / (1 + RESAMPLES)
        dof = len(t) - design.shape[1] - 2
        p_gauss = min(1., len(freq) * f_dist.sf(peak * dof / (2 * max(1 - peak, 1e-14)), 2, dof))
        seed = freq[idx[0]]
        p, k, e, omega, tp, coef, res, bases, success = fit_orbit(t, y, err, design, seed, span)
        frac = ((t - tp) / p) % 1
        phase_gap = float(np.diff(np.sort(np.r_[frac, frac.min() + 1])).max())
        sse0, sse = np.sum((residual / err)**2), res @ res
        delta_bic = len(t) * np.log(max(sse0, 1e-30) / max(sse, 1e-30)) - 5 * np.log(len(t))
        # Subtract previously fitted companions (at their fixed orbital shapes).
        y_current = y - design[:, base.shape[1]:] @ coef[base.shape[1]:-2]
        activity_info, flags = activity(n, y_current, 1 / p, span)
        competing, at_window = aliases(freq, powers, windows, span, 1 / p)
        alias_fits = []
        if fap <= .01 and delta_bic >= 10:
            for rival in competing[:3]:
                ap, ak, ae, _, _, _, ares, _, ok = fit_orbit(t, y, err, design, 1 / rival, span)
                alias_fits.append(dict(period_d=ap, k_ms=ak, ecc=ae, converged=bool(ok),
                                       bic_minus_primary=len(t) * np.log(max(ares @ ares, 1e-30) / max(sse, 1e-30))))
        shape = bases @ np.array([np.cos(omega), -np.sin(omega)])
        split, split_flags = amplitude_checks(t, y_current, err, era, shape)
        flags.extend(split_flags)
        if competing or at_window:
            flags.append("window_alias")
        if phase_gap > .35 or span / p < 1.5:
            flags.append("incomplete_coverage")
        if not success or e > .82:
            flags.append("orbit_poorly_constrained")
        if n_blocks < 8:
            flags.append("too_few_noise_blocks")
        if fap > .01:
            flags.append("bootstrap_fap")
        survey_bound = min(1., 3 * star_count * p_gauss)
        if survey_bound >= .05:
            flags.append("survey_trials")
        if delta_bic < 10:
            flags.append("weak_orbit")
        reject = ("bootstrap_fap" in flags or "weak_orbit" in flags or
                  any(v.startswith("activity_") or v.startswith("rotation_proxy_") for v in flags))
        status = "rejected" if reject else ("inconclusive" if flags else "pass")
        signals.append(dict(star=star, signal_index=seq, period_d=p, k_ms=k, ecc=e,
                            omega_rad=omega, tperi_bjd=tp, grid_period_d=1 / seed,
                            fap=fap, bootstrap_exceedances=exc, n_bootstrap=RESAMPLES,
                            fap_gaussian_grid=p_gauss, gaussian_trials_screen=survey_bound,
                            n_blocks=n_blocks, peak_power=peak, delta_bic=delta_bic,
                            phase_gap=phase_gap, offset_post_minus_pre_ms=coef[2] if base.shape[1] == 3 else "",
                            window_peaks_d=json.dumps((1 / windows).round(6).tolist()),
                            alias_periods_d=json.dumps(competing), alias_fits_json=json.dumps(alias_fits),
                            activity_json=json.dumps(activity_info),
                            amplitude_checks_json=json.dumps(split), rv_status=status,
                            flags="|".join(flags)))
        # Model significant components, even if associated with activity or a known binary.
        # Keep both quadratures free when fitting additional components.
        if fap > .01 or delta_bic < 10 or not success or n_blocks < 8:
            break
        extended = np.column_stack([design, bases])
        if np.linalg.cond(extended / err[:, None]) > 1e8:
            break
        design = extended
        if len(t) <= design.shape[1] + 10:
            break
    summary["n_scans"] = len(signals)
    coef, std = trend(t, y, err, design, blocks)
    summary.update(slope_ms_yr=coef[1], slope_err_ms_yr=std[1])
    # Calibrate the residual detection threshold on the noise used for the injections.
    final_grid = Grid(t, err, design, freq)
    residual = final_grid.project(y)[0] * err
    if design.shape[1] != grid.q.shape[1]:
        signs = rng.integers(0, 2, (RESAMPLES, n_blocks)) * 2 - 1
        boot, _, _ = final_grid.scan(residual * signs[:, blocks])
    else:
        boot = high[1:]
    threshold = np.sort(boot)[-2]
    return summary, signals, inject(star, final_grid, residual, threshold, rng)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", required=True)
    parser.add_argument("--targets")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--shard", type=int, default=0)
    args = parser.parse_args()
    if os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_EVENT_NAME") != "workflow_dispatch":
        raise SystemExit("This search runs in the manually dispatched GitHub workflow.")
    if not (0 <= args.shard < args.shards):
        raise ValueError("invalid shard index")
    df = pd.read_csv(args.epochs, dtype={"star": str}, low_memory=False)
    groups = df.groupby("star", sort=False)
    stars = set(groups.groups)
    if args.targets:
        stars.update(pd.read_csv(args.targets, usecols=["star"], dtype=str).star.dropna())
    total = len(stars)
    args.out.mkdir(parents=True, exist_ok=True)
    failures = 0
    with ExitStack() as files:
        writers = []
        for name, cols in [("stars", STAR_COLUMNS), ("signals", SIGNAL_COLUMNS), ("injections", INJECTION_COLUMNS)]:
            file = files.enter_context((args.out / (name + ".csv")).open("w", newline="", encoding="utf-8"))
            writer = csv.DictWriter(file, fieldnames=cols)
            writer.writeheader()
            writers.append(writer)
        for star in sorted(stars)[args.shard::args.shards]:
            frame = groups.get_group(star) if star in groups.groups else df.iloc[:0]
            try:
                summary, signals, injections = search(star, frame, total)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
                failures += 1
                summary = dict(star=star, n_epochs=len(frame), eligibility="error", error=str(exc))
                signals, injections = [], []
            writers[0].writerow(summary)
            writers[1].writerows(signals)
            writers[2].writerows(injections)
    print(f"Processed shard {args.shard}/{args.shards}, per-target errors={failures}")
    if failures:
        raise SystemExit("Inspect and resolve per-target errors before claiming complete coverage.")


if __name__ == "__main__":
    main()
