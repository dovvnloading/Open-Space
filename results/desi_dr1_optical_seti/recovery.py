"""Injection/recovery trials for one star and arm. Keep the full trial ledger."""
from dataclasses import asdict, dataclass

import numpy as np
from scipy.stats import beta

from laser_search import Settings, inject_line, scan_target


@dataclass(frozen=True)
class Trial:
    seed: int
    index: int
    rest_A: float
    flux17: float
    eligible: bool
    recovered: bool
    reason: str


def _matches(candidate, wavelength):
    # The injected wavelength is exact. Agreement must hold in every exposure.
    return all(abs(h.rest_A - wavelength) <= min(0.5*h.rest_fwhm_A, 3*h.rest_error_A)
               for h in candidate.hits)


def measure_completeness(exposures, expected_expids, rest_range, fluxes, trials, seed,
                         settings=Settings(), refit=None, veto=None):
    """Return trial records for random uniform vacuum rest wavelengths.

    Same draws are paired across fluxes. Masks/gaps count as misses in the full
    interval. Fixed-template runs include the exposure-local affine fit. Supply
    refit(injected_exposures) to propagate line absorption into stellar fitting.
    veto(exposures,candidate) returns None for a pass or a reason string.
    Omitting either callback limits inference to its corresponding pipeline stage.
    Detector shape/persistence is always rerun.
    """
    exposures, expected_expids = tuple(exposures), tuple(expected_expids)
    baseline = scan_target(exposures, expected_expids, settings)
    lo, hi = rest_range
    if not np.all(np.isfinite([lo, hi])) or not (0 < lo < hi) or trials < 1:
        raise ValueError("need a positive increasing rest range and positive trial count")
    fluxes = tuple(float(value) for value in fluxes)
    if not fluxes or len(set(fluxes)) != len(fluxes) or not all(np.isfinite(f) and f > 0 for f in fluxes):
        raise ValueError("flux grid must have distinct positive finite integrated fluxes")
    wavelengths = np.random.default_rng(seed).uniform(lo, hi, trials)
    records = []
    for i, wave in enumerate(wavelengths):
        eligible = any(a <= wave <= b for a, b in baseline.coverage)
        # Never credit a line already seen before injection, even if astrophysical
        # vetting would subsequently reject it. Report confusion explicitly.
        confused = any(_matches(c, wave) for c in baseline.candidates)
        for f in fluxes:
            if not eligible:
                reason = "outside_joint_coverage"
            elif confused:
                reason = "baseline_coincidence"
            else:
                injected = tuple(inject_line(e, wave, f) for e in exposures)
                if refit is not None:
                    injected = tuple(refit(injected))
                search = scan_target(injected, expected_expids, settings)
                near = [c for c in search.candidates if _matches(c, wave)]
                reason = "missed_detection"
                for candidate in near:
                    rejection = veto(injected, candidate) if veto is not None else None
                    if rejection is None:
                        reason = "recovered"
                        break
                    reason = "veto_" + str(rejection)
            records.append(asdict(Trial(int(seed), i, float(wave), f, bool(eligible), reason == "recovered", reason)))
    return records


def summarize_recovery(trials, fraction=0.95, confidence=0.95, eligible_only=False):
    """Summarize the ledger and obtain a conservative grid-point flux threshold.

    Use an exact one-sided binomial lower bound with Bonferroni correction over
    all tested levels. Every tested level from the reported threshold up must
    clear the requested recovery fraction. No extrapolation past tested levels.
    Do not pool correlated repeats of the same injection draws or different stars.
    """
    if not (0 < fraction < 1 and 0 < confidence < 1):
        raise ValueError("fraction and confidence must be between zero and one")
    fluxes = sorted(set(t["flux17"] for t in trials))
    if not fluxes:
        raise ValueError("empty trial ledger")
    alpha = (1. - confidence) / len(fluxes)
    summary = []
    for flux in fluxes:
        rows = [row for row in trials if row["flux17"] == flux and (row["eligible"] or not eligible_only)]
        keys = [(r["seed"],r["index"]) for r in rows]
        if len(keys) != len(set(keys)):
            raise ValueError("duplicate trials at the same flux")
        n = len(rows)
        k = sum(r["recovered"] for r in rows)
        lower = float(beta.ppf(alpha, k, n-k+1)) if k else 0.
        summary.append(dict(flux17=flux, trials=n, recovered=k, fraction=k/n if n else None, lower=lower))
    limit = None
    for row in reversed(summary):
        if row["lower"] < fraction:
            break
        limit = row["flux17"]
    return summary, limit
