"""Sampling, model-dependence and interval checks for OGLE-SMC-CEP-2332."""

import csv
import json
import numpy as np
from scipy.stats import spearmanr
from screen import ROOT, catalog, read_curves, basis, estimate
from modulation import segments, model, design, score
from waveform import delve

STAR = "OGLE-SMC-CEP-2332"


def interval_checks(curves, p, q, rng):
    result = {}
    for key, left in [("weak", 5510), ("strong", 6080)]:
        # Same windows as 190-d screen, selected by measured minimum/maximum of the later phase.
        data = curves[STAR, 4, "I"]
        a = data[np.floor(data[:, 0] / 190) == np.floor(left / 190)]
        n = np.floor(a[:, 0]).astype(int)
        nights = np.unique(n)
        draws = []
        freqs = np.linspace(1 / q - 0.003, 1 / q + 0.003, 151)
        # Refit both overtone phases and harmonics independently in each resample.
        groups = [a[n == k] for k in nights]
        for i in range(800):
            b = estimate(
                np.concatenate(
                    [groups[k] for k in rng.integers(len(nights), size=len(nights))]
                ),
                p,
                q,
            )
            if b is not None:
                draws.append(b)
        # Test loss of phase coherence across one interval by allowing the second frequency to float.
        search = []
        for h in freqs:
            X = basis(a[:, 0], p, 1 / h)
            w = 1 / np.hypot(a[:, 2], 0.009)
            coeff = np.linalg.lstsq(X * w[:, None], a[:, 1] * w, rcond=None)[0]
            err = np.sum(((a[:, 1] - X @ coeff) * w) ** 2)
            search.append((err, h, float(np.hypot(*coeff[3:5]))))
        err, h, refit = min(search)
        result[key] = dict(
            n=len(a),
            nights=len(nights),
            used_draws=len(draws),
            avg_time=float(a[:, 0].mean()),
            date_range=[float(a[:, 0].min()), float(a[:, 0].max())],
            A1_A2=estimate(a, p, q),
            A2_bootstrap=np.percentile(np.array(draws)[:, 1], [2.5, 50, 97.5]).tolist(),
            free_frequency_A2=refit,
            best_freq=h,
        )
    return result


def phase_offset_check(seq_early, seq_late, number=2000):
    # For fixed early phases, compare a rotation against real held-out values. Tracks phase specificity.
    t, y = seq_late[:, 0], seq_late[:, 2]
    P, b = model(seq_early)
    mse = float(np.mean((y - design(t, P) @ b) ** 2))
    phases = np.arange(number) / number
    errors = [
        float(np.mean((y - design(t + P * phase, P) @ b) ** 2)) for phase in phases
    ]
    # Describes agreement for a grid of offsets; no false-alarm claim after ranking and target selection.
    return dict(
        period=P,
        model_beta=b.tolist(),
        later_measured_mse=mse,
        fraction_phase_offsets_as_good=float(np.mean(np.array(errors) <= mse)),
        holdout=score(seq_late, P, b),
    )


def decam_check():
    data, bands = delve("2332")
    p, q = catalog()[STAR][1:]
    # The isolated, dense 2020 visit consists of four filters. Require a full-rank harmonic design.
    selected = (data[:, 0] > 9140) & (data[:, 0] < 9170)
    report = {}
    for band in "griz":
        a = data[selected & (bands == band)]
        t, y, e = a.T
        phase = 2 * np.pi * (t - np.median(t))
        X = np.column_stack(
            [np.ones(len(t))]
            + [
                fun(phase * freq)
                for freq in [1 / p, 2 / p, 1 / q]
                for fun in [np.cos, np.sin]
            ]
        )

        def fit(mask):
            xx = X[mask]
            yy = y[mask]
            c = np.linalg.lstsq(xx, yy, rcond=None)[0]
            return (
                float(np.hypot(*c[1:3])),
                float(np.hypot(*c[-2:])),
                float(np.linalg.cond(xx)),
            )

        n = np.floor(t)
        jack = [fit(n != d)[1] for d in np.unique(n)]
        A1, A2, cond = fit(np.ones(len(t), bool))
        report[band] = dict(
            n=len(t),
            nights=len(set(n)),
            t_min=float(t.min()),
            t_max=float(t.max()),
            A1=A1,
            A2=A2,
            condition=cond,
            leave_night_out_A2_range=[min(jack), max(jack)],
        )
    return report


def injection(data, p, q, rng, ntrial=350):
    # Constant-amplitude artificial pulsations, sampled on the real cadence. The residuals include
    # actual modulation as extra noise. Resample whole nights to preserve within-night correlations.
    all_windows = []
    for a in data:
        bins = np.floor(a[:, 0] / 190)
        for k in np.unique(bins):
            win = a[bins == k]
            if (
                len(win) >= 24
                and np.ptp(win[:, 0]) >= 28
                and estimate(win, p, q) is not None
            ):
                all_windows.append(win)
    norms = []
    amp0 = 0.03
    for win in all_windows:
        t, y, e = win.T
        X = basis(t, p, q)
        w = 1 / np.hypot(e, 0.007)
        c = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)[0]
        resid = y - X @ c
        c[3:5] *= amp0 / max(np.hypot(*c[3:5]), 1e-8)
        model_y = X @ c
        _, which = np.unique(np.floor(t), return_inverse=True)
        sims = []
        for _ in range(ntrial):
            sign = rng.choice([-1.0, 1.0], size=which.max() + 1)
            z = estimate(np.column_stack([t, model_y + sign[which] * resid, e]), p, q)
            sims.append(z[1] if z is not None else np.nan)
        norms.append(sims)
    peak_to_peak = np.nanmax(norms, axis=0) - np.nanmin(norms, axis=0)
    return dict(
        injected_A2=amp0,
        ntrial=ntrial,
        nwindows=len(all_windows),
        p2p_percentiles=np.percentile(peak_to_peak, [50, 95, 99, 100]).tolist(),
    )


def main():
    rng = np.random.default_rng(813507)
    p, q = catalog()[STAR][1:]
    curves = dict(read_curves([STAR]))
    early = curves[STAR, 3, "I"]
    late = curves[STAR, 4, "I"]
    a = segments(early, p, q, 190, 0)
    b = segments(late, p, q, 190, 0)
    s = np.vstack([a, b])
    P, beta = model(s)
    leave = [model(np.delete(s, i, axis=0))[0] for i in range(len(s))]
    phase = np.arange(0, P, 0.1)
    mag = design(phase + 3500, P) @ beta
    out = dict(
        early_range=early[:, 0][
            [np.argmin(early[:, 0]), np.argmax(early[:, 0])]
        ].tolist(),
        late_range=late[:, 0][[np.argmin(late[:, 0]), np.argmax(late[:, 0])]].tolist(),
        n_early=len(early),
        n_late=len(late),
        n_bins_early=len(a),
        n_bins_late=len(b),
        combined_period=P,
        leave_bin_out_period_range=[min(leave), max(leave)],
        phase_max_HJD_minus_2450000=3500 + float(phase[np.argmax(mag)]),
        phase_min_HJD_minus_2450000=3500 + float(phase[np.argmin(mag)]),
        spearman=float(spearmanr(s[:, 1], s[:, 2]).statistic),
        split=phase_offset_check(a, b),
        intervals=interval_checks(curves, p, q, rng),
        cadence=injection([early, late], p, q, rng),
        decam_2020=decam_check(),
    )
    (ROOT / "validation_2332.json").write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
