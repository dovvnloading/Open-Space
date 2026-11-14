"""Cross-predict whole seasons of eclipse times and propagate the Gaia orbit covariance."""

import argparse
import ast
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

ROOT = Path(__file__).resolve().parent
FIELDS = ["c_thiele_innes", "h_thiele_innes", "eccentricity", "period", "t_periastron"]
INDICES = [9, 10, 12, 13, 14]  # AstroSpectroSB1 full 15-parameter order
REF_JD = Time(2016.0, format="jyear", scale="tcb").jd
C_AU_S = 499.00478384


def covariance(row):
    if row.bit_index != 65535:
        raise ValueError("Expected the full AstroSpectroSB1 solution")
    packed = np.asarray(ast.literal_eval(row.corr_vec))
    corr = np.eye(15)
    # Column-major upper triangle from the Gaia DR3 data model.
    for j in range(1, 15):
        for i in range(j):
            corr[i, j] = corr[j, i] = packed[j * (j - 1) // 2 + i]
    err = np.array([row[key + "_error"] for key in FIELDS])
    return corr[np.ix_(INDICES, INDICES)] * np.outer(err, err)


def delays(t, orbit, shift=0):
    c, h, ecc, period, periastron = orbit
    mean_anomaly = 2 * np.pi * (((t - periastron) / period + shift) % 1)
    E = mean_anomaly.copy()
    for _ in range(10):
        E -= (E - ecc * np.sin(E) - mean_anomaly) / (1 - ecc * np.cos(E))
    return C_AU_S * (c * (np.cos(E) - ecc) + h * np.sqrt(1 - ecc**2) * np.sin(E))


def design_matrix(t, degree, extras):
    u = (t - 2750) / 1000
    return np.column_stack([u**k for k in range(degree + 1)] + list(extras))


def cv(t, y, groups, degree, extras, mode="absolute"):
    """All model coefficients fitted on other visits. Last two sectors form one visit."""
    x = design_matrix(t, degree, extras)
    errors = []
    for group in np.unique(groups):
        held = groups == group
        train = ~held
        xtrain, ytrain = x[train].copy(), y[train].copy()
        xtest, ytest = x[held].copy(), y[held].copy()
        if mode in ("first-anchored", "centered-shape"):
            # Profile out a constant on each TRAINING visit with its full training data.
            for g in np.unique(groups[train]):
                select = groups[train] == g
                xtrain[select] -= xtrain[select].mean(0)
                ytrain[select] -= ytrain[select].mean()
            xtrain = xtrain[:, 1:]
            if mode == "first-anchored":
                # Pay one held-out clock measurement to set its offset. Score the remainder.
                ytest = ytest[1:] - ytest[0]
                xtest = xtest[1:, 1:] - xtest[0, 1:]
            else:
                # Post-hoc shape diagnostic only: discards all held-out zero-point information.
                ytest -= ytest.mean()
                xtest = xtest[:, 1:] - xtest[:, 1:].mean(0)
        beta = np.linalg.lstsq(xtrain, ytrain, rcond=None)[0]
        residual = ytest - xtest @ beta
        errors.append(np.sqrt(np.mean(residual**2)))
    return np.array(errors)


def equal_visit(rms):
    return np.sqrt(np.mean(np.square(rms), axis=-1))


def sinusoid(t, period):
    return [np.sin(2 * np.pi * t / period), np.cos(2 * np.pi * t / period)]


def main(source, dest, draws):
    dest.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(source / "outputs/timings.csv").dropna(subset=["both"])
    orbit = pd.read_csv(source / "data/gaia-orbit.csv").iloc[0]
    times = table.time.to_numpy()
    y = table.both.to_numpy()
    sectors = table.sector.to_numpy()
    visits = np.minimum(sectors, 103)
    tcb = Time(times + 2457000, format="jd", scale="tdb").tcb.jd - REF_JD
    params = orbit[FIELDS].astype(float).to_numpy()
    cov = covariance(orbit)
    z = delays(tcb, params)
    add = {
        "none": [],
        "Gaia": [z],
        "annual": sinusoid(times, 365.25),
        "semiannual": sinusoid(times, 365.25 / 2),
        "free_sine_at_Gaia_period": sinusoid(times, orbit.period),
    }
    folds = []
    for clock in ["both", "primary", "secondary", "primary_minus_secondary"]:
        yy = (
            y
            if clock == "both"
            else table.primary.to_numpy() - table.secondary.to_numpy()
            if clock == "primary_minus_secondary"
            else table[clock].to_numpy()
        )
        for degree in [1, 2]:
            for name, extra in add.items():
                if clock != "both" and name not in ("none", "Gaia"):
                    continue
                for mode in ("absolute", "first-anchored", "centered-shape"):
                    scores = cv(times, yy, visits, degree, extra, mode)
                    folds.append(
                        dict(
                            clock=clock,
                            degree=degree,
                            term=name,
                            mode=mode,
                            equal_visit_rms=equal_visit(scores),
                            **dict(zip(("S12", "S39", "S66", "S103_104"), scores)),
                        )
                    )
    pd.DataFrame(folds).to_csv(
        dest / "season-predictions.csv", index=False, float_format="%.8f"
    )
    controls = []
    for phase in np.arange(100) / 100:
        zp = delays(tcb, params, phase)
        for degree in [1, 2]:
            for mode in ["absolute", "first-anchored"]:
                controls.append(
                    (
                        phase,
                        degree,
                        mode,
                        equal_visit(cv(times, y, visits, degree, [zp], mode)),
                    )
                )
    ctr = pd.DataFrame(
        controls, columns=["phase_shift", "degree", "mode", "equal_visit_rms"]
    )
    ctr.to_csv(dest / "phase-controls.csv", index=False, float_format="%.8f")

    # Preserve every cadence within its visit. Permute only which of the four Gaia-predicted
    # local delay slopes is assigned to each observed visit; this is a diagnostic, not a FAP.
    slope_rows = []
    for g in np.unique(visits):
        sel = visits == g
        observed = {
            k: np.polyfit(times[sel] - times[sel].mean(), table.loc[sel, k], 1)[0]
            for k in ("both", "primary", "secondary")
        }
        slope_rows.append(
            dict(
                visit=int(g),
                predicted=np.polyfit(times[sel] - times[sel].mean(), z[sel], 1)[0],
                **observed,
            )
        )
    slopes = pd.DataFrame(slope_rows)
    slopes.to_csv(dest / "visit-slopes.csv", index=False, float_format="%.8f")
    inds = range(len(slopes))
    permutations = []
    for order in itertools.permutations(inds):
        xx = np.column_stack(
            [np.ones(len(inds)), slopes.predicted.to_numpy()[list(order)]]
        )
        yy = slopes.both.to_numpy()
        b = np.linalg.lstsq(xx, yy, rcond=None)[0]
        permutations.append(
            {
                "assignment": "/".join(str(slopes.visit.iloc[i]) for i in order),
                "scale": b[-1],
                "rms_seconds_per_day": np.mean((yy - xx @ b) ** 2) ** 0.5,
            }
        )
    pd.DataFrame(permutations).to_csv(
        dest / "slope-permutations.csv", index=False, float_format="%.8f"
    )
    periods = sorted(
        set(np.arange(120, 601, 2).astype(float)) | {params[3], 365.25, 365.25 / 2}
    )
    period_rows = []
    # Keep periastron epoch, orientation and eccentricity fixed while varying period.
    # This particular grid is only a specificity diagnostic of the sparse cadence.
    for trial in periods:
        varied = params.copy()
        varied[3] = trial
        alt = delays(tcb, varied)
        for degree in (1, 2):
            period_rows.append(
                {
                    "trial_period_days": trial,
                    "degree": degree,
                    "equal_visit_rms_s": equal_visit(
                        cv(times, y, visits, degree, [alt])
                    ),
                }
            )
    pd.DataFrame(period_rows).to_csv(
        dest / "period-controls.csv", index=False, float_format="%.8f"
    )

    rng = np.random.default_rng(96100077)
    cloud = rng.multivariate_normal(params, cov, size=draws)
    signals = delays(tcb[None, :], cloud.T[:, :, None])
    quantiles = [0.025, 0.5, 0.975]
    mc = {}
    for degree in (1, 2):
        for mode in ("absolute", "first-anchored"):
            scores = np.array(
                [cv(times, y, visits, degree, [zz], mode) for zz in signals]
            )
            null = equal_visit(cv(times, y, visits, degree, [], mode))
            mc[f"degree_{degree}_{mode}"] = {
                "equal_visit_rms_quantiles_s": np.quantile(
                    equal_visit(scores), quantiles
                ).tolist(),
                "fraction_lower_rms_than_same_baseline": float(
                    np.mean(equal_visit(scores) < null)
                ),
                "forward_2026_rms_quantiles_s": np.quantile(
                    scores[:, -1], quantiles
                ).tolist(),
            }
    rows = []
    for n in range(13, 18):
        epoch = REF_JD + cloud[:, 4] + n * cloud[:, 3]
        interval = Time(
            np.quantile(epoch, [0.025, 0.16, 0.5, 0.84, 0.975]),
            scale="tcb",
            format="jd",
        ).utc.isot
        rows.append((n, *interval))
    pd.DataFrame(
        rows, columns=["cycle", "utc_p025", "utc_p16", "utc_p50", "utc_p84", "utc_p975"]
    ).to_csv(dest / "periastra.csv", index=False)
    (dest / "covariance.json").write_text(
        json.dumps(
            {
                "seed": 96100077,
                "draws": draws,
                "quantiles": quantiles,
                "fields": FIELDS,
                "parameter_covariance": cov.tolist(),
                "season_forecasts": mc,
            },
            indent=2,
        )
        + "\n"
    )
    print(
        pd.DataFrame(folds)
        .query("clock=='both' and mode=='absolute'")[
            ["degree", "term", "equal_visit_rms", "S103_104"]
        ]
        .round(1)
        .to_string(index=False)
    )
    print(mc)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path, help="original results/tic-96100077 directory")
    p.add_argument("--output", type=Path, default=ROOT / "outputs")
    p.add_argument("--draws", type=int, default=2000)
    args = p.parse_args()
    main(args.source, args.output, args.draws)
