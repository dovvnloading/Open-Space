"""Diagnostic harmonic comparisons on the observed nights; no candidate cuts.

For measured f and a reference fundamental f0, the signed separation in
Rayleigh widths is (f - n*f0)*baseline. Window amplitude is the absolute
weighted mean of exp(2*pi*i*f*t), normalized to one at zero frequency. The
canonical correlation is the maximum over phases of the two sampled sinusoid
spaces after weighting and removing the common offset, drift and fibre terms.
A high correlation describes the cadence's ability to distinguish periods.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from runtime import require_survey_workflow
from search import nightly, nuisance, orbit_basis

YEAR = 365.25
FIELDS = """star signal_index period_d k_ms fap delta_bic rv_status reference_type
reference_signal_index reference_period_d reference_fap reference_delta_bic harmonic
harmonic_period_d signed_rayleigh abs_rayleigh within_one_rayleigh n_nights baseline_d
window_candidate window_harmonic window_difference window_sum canonical_correlation
reference_orbit_projection_k_ms""".split()
SUMMARY = """star signal_index period_d annual_harmonic annual_abs_rayleigh annual_window_amplitude
annual_canonical_correlation other_significant_index other_fundamental_d other_harmonic
other_abs_rayleigh other_canonical_correlation other_orbit_projection_k_ms""".split()


def orthogonal_basis(matrix):
    u, scales, _ = np.linalg.svd(matrix, full_matrices=False)
    keep = scales > np.finfo(float).eps * max(matrix.shape) * scales[0]
    return u[:, keep]


def run(work, out):
    require_survey_workflow()
    fits = pd.read_csv(out / "signals.csv", dtype={"star": str})
    selected = fits[(fits.fap <= .01) & (fits.delta_bic >= 10)]
    epochs = pd.read_csv(work / "epochs.csv.gz", dtype={"star": str})
    grouped = epochs[epochs.star.isin(selected.star)].groupby("star", sort=False)
    comparisons, summary = [], []
    for star, tested in selected.groupby("star", sort=True):
        obs = nightly(grouped.get_group(star))
        err = obs.err.to_numpy()
        t = obs.bjd.to_numpy() - obs.bjd.min()
        duration = np.ptp(t)
        weights = err**-2 / np.sum(err**-2)
        q = orthogonal_basis(nuisance(t, obs.era.to_numpy()) / err[:, None])

        def window(f):
            return float(abs(np.sum(weights * np.exp(2j*np.pi*f*t))))

        def design(freq):
            phase = 2*np.pi*freq*t
            a = np.column_stack([np.sin(phase), np.cos(phase)]) / err[:, None]
            return a - q @ (q.T @ a)

        for sig in tested.itertuples(index=False):
            freq = 1 / sig.period_d
            fit_design = design(freq)
            basis = orthogonal_basis(fit_design)
            observed_window = window(freq)
            refs = [("year", None, YEAR, np.nan, np.nan, np.nan)]
            others = fits[fits.star.eq(star) & fits.signal_index.ne(sig.signal_index)]
            for row in others.itertuples(index=False):
                # How much of the other fitted Keplerian projects onto a sine
                # at this candidate period? This includes window leakage.
                m0 = 2*np.pi*(obs.bjd.min() - row.tperi_bjd) / row.period_d
                orbit = orbit_basis(t, 1/row.period_d, row.ecc, m0)
                model = row.k_ms * (orbit @ [np.cos(row.omega_rad), -np.sin(row.omega_rad)]) / err
                beta = np.linalg.lstsq(fit_design, model - q @ (q.T @ model), rcond=None)[0]
                refs.append(("fitted orbit", row.signal_index, row.period_d, row.fap,
                             row.delta_bic, float(np.linalg.norm(beta))))
            closest = {}
            for kind, index, period, ref_fap, ref_bic, leakage in refs:
                for n in range(1, 7):
                    fref = n/period
                    separation = (freq - fref)*duration
                    other_basis = orthogonal_basis(design(fref))
                    cross = basis.T @ other_basis
                    cc = min(1., np.linalg.norm(cross, ord=2)) if cross.size else np.nan
                    row = dict(star=star, signal_index=sig.signal_index, period_d=sig.period_d,
                               k_ms=sig.k_ms, fap=sig.fap, delta_bic=sig.delta_bic, rv_status=sig.rv_status,
                               reference_type=kind, reference_signal_index=index, reference_period_d=period,
                               reference_fap=ref_fap, reference_delta_bic=ref_bic, harmonic=n,
                               harmonic_period_d=1/fref, signed_rayleigh=separation, abs_rayleigh=abs(separation),
                               within_one_rayleigh=abs(separation)<=1, n_nights=len(t), baseline_d=duration,
                               window_candidate=observed_window, window_harmonic=window(fref),
                               window_difference=window(freq-fref), window_sum=window(freq+fref),
                               canonical_correlation=cc, reference_orbit_projection_k_ms=leakage)
                    comparisons.append(row)
                    if kind == "year" or (ref_fap <= .01 and ref_bic >= 10):
                        if kind not in closest or abs(separation) < closest[kind]["abs_rayleigh"]:
                            closest[kind] = row
            year = closest["year"]
            best = dict(star=star, signal_index=sig.signal_index, period_d=sig.period_d,
                        annual_harmonic=year["harmonic"], annual_abs_rayleigh=year["abs_rayleigh"],
                        annual_window_amplitude=year["window_harmonic"],
                        annual_canonical_correlation=year["canonical_correlation"])
            if "fitted orbit" in closest:
                other = closest["fitted orbit"]
                best.update(other_significant_index=other["reference_signal_index"],
                            other_fundamental_d=other["reference_period_d"], other_harmonic=other["harmonic"],
                            other_abs_rayleigh=other["abs_rayleigh"], other_canonical_correlation=other["canonical_correlation"],
                            other_orbit_projection_k_ms=other["reference_orbit_projection_k_ms"])
            summary.append(best)
    pd.DataFrame(comparisons, columns=FIELDS).to_csv(out / "alias-review.csv", index=False)
    pd.DataFrame(summary, columns=SUMMARY).to_csv(out / "alias-review-summary.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    run(args.work, args.out)
