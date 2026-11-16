"""Reproducible independent synthetic controls for the O4b H1/L1 burst scan."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
import scipy
from scipy import fft, signal, stats

import core

FS = 2048
SECONDS = 160  # 16 s of extra guard on either end of a 128 s experiment
ON_START, ON_STOP = 48, 112
SHIFTS = tuple(range(-32, 33, 4))
SEED = 405116
N_REALIZATIONS = 192
CENTER = 64.31
DELAY = 8 / FS
H_RSS = 2e-21


def colored_noise(rng, detector):
    # The target ASD is physical strain/sqrt(Hz). The records are independent
    # across sites; each has a different low-frequency knee and high-frequency rise.
    level, low, high = {"H1": (3e-23, 55, 330), "L1": (2.5e-23, 46, 410)}[detector]
    n = FS * SECONDS
    freq = fft.rfftfreq(n, 1 / FS)
    asd = level * np.sqrt((1 + (low / np.maximum(freq, 12)) ** 4) * (1 + (freq / high) ** 2))
    return fft.irfft(fft.rfft(rng.normal(size=n)) * np.sqrt(FS / 2) * asd, n=n)


def run_search(h, l, shifts):
    valid = np.ones(SECONDS, dtype=bool)
    a, b = core.prepare(h, FS, valid), core.prepare(l, FS, valid)
    # The primary-PSD training can use independent flanks, also at zero lag.
    a.live[:ON_START] = False
    a.live[ON_STOP:] = False
    result = core.scan(a, b, shifts)
    if any(r.live_seconds != ON_STOP - ON_START for r in result.values()):
        raise ValueError("Unequal simulated exposure across slides")
    return result


def max_score(slide):
    return max((t.score for t in slide.triggers), default=0.0)


def make_test_wave(rng):
    # A Gaussian-enveloped broadband random transient, independent of the
    # analysis statistic and its sine-Gaussian sensitivity injections.
    x = np.arange(FS * SECONDS) / FS - CENTER
    band_filter = signal.butter(5, [96, 224], btype="bandpass", fs=FS, output="sos")
    band_noise = signal.sosfiltfilt(band_filter, rng.normal(size=x.size))
    pulse = band_noise * np.exp(-0.5 * (x / 0.035) ** 2)
    return pulse * H_RSS / np.sqrt(np.sum(pulse**2) / FS)


def injection_tests(seed):
    hseed, lseed, burstseed = seed.spawn(3)
    h = colored_noise(np.random.default_rng(hseed), "H1")
    l = colored_noise(np.random.default_rng(lseed), "L1")
    burst = make_test_wave(np.random.default_rng(burstseed))
    shifted = np.zeros_like(burst)
    shifted[8:] = -0.8 * burst[:-8]
    cases = {"clean_noise": (h, l), "hanford_only": (h + burst, l),
             "coherent": (h + burst, l + shifted)}
    results = {}
    for name, (a, b) in cases.items():
        slide = run_search(a, b, (0,))[0]
        peak = max(slide.triggers, key=lambda t: t.score, default=None)
        nearby = [t for t in slide.triggers if abs(t.time - CENTER) <= 0.5]
        results[name] = {"loudest": asdict(peak) if peak else None,
                         "target_score": max((t.score for t in nearby), default=0.0)}
    peak = results["coherent"]["loudest"]
    return {"center_s": CENTER, "L_minus_H_s": DELAY, "L_over_H": -0.8,
            "H_hrss_strain_sqrt_s": H_RSS, "recovery_time_tolerance_s": 0.5,
            "recovery_lag_tolerance_s": 0.0015, "cases": results,
            "recovered_location_lag_sign": (
                peak is not None and abs(peak["time"] - CENTER) <= 0.5 and
                abs(peak["lag_seconds"] - DELAY) <= 0.0015 and peak["sign"] == -1)}


def main(output):
    child_seeds = np.random.SeedSequence(SEED).spawn(N_REALIZATIONS + 1)
    zero_idx = SHIFTS.index(0)
    all_maxima = []
    ntiles, gated = 0, 0
    for seed in child_seeds[:-1]:
        hs, ls = seed.spawn(2)
        h = colored_noise(np.random.default_rng(hs), "H1")
        l = colored_noise(np.random.default_rng(ls), "L1")
        scans = run_search(h, l, SHIFTS)
        all_maxima.append([max_score(scans[s]) for s in SHIFTS])
        ntiles += sum(r.nominal_tiles for r in scans.values())
        gated += sum(r.gated_tiles for r in scans.values())
    maxima = np.array(all_maxima)
    on = maxima[:, zero_idx]
    off = np.delete(maxima, zero_idx, axis=1)
    # An on/off rank is descriptive: off-time comparisons reuse the same two
    # records, and alignment with each local PSD training half can matter.
    ranks = 1 + np.sum(off >= on[:, None], axis=1)
    # The first 96 zero-lag maxima set a threshold before independent validation.
    # For IID continuous maxima the held-out exceedance count is beta-binomial.
    n = N_REALIZATIONS // 2
    k = int(np.ceil(0.95 * (n + 1)))  # 1-based calibration order statistic
    threshold = float(np.sort(on[:n])[k - 1])
    validation = int(np.sum(on[n:] > threshold))
    prediction = stats.betabinom.interval(0.95, n, n + 1 - k, k)
    report = {
        "scope": "synthetic only; no observing time or search-background exposure",
        "seed": SEED, "numpy_version": np.__version__, "scipy_version": scipy.__version__,
        "sample_rate": FS, "generated_seconds_per_detector": SECONDS,
        "search_seconds_per_realization": ON_STOP - ON_START, "realizations": N_REALIZATIONS,
        "shifts_s": SHIFTS, "foreground_seconds": N_REALIZATIONS * (ON_STOP - ON_START),
        "offtime_reused_seconds": N_REALIZATIONS * (ON_STOP - ON_START) * (len(SHIFTS) - 1),
        "unshifted_and_shifted_nominal_tiles": ntiles, "double_gate_tiles": gated,
        "global_rank_histogram_1_is_loudest": np.bincount(ranks, minlength=len(SHIFTS) + 1)[1:].tolist(),
        "foreground_maximum_median": float(np.median(on)),
        "offtime_maximum_median": float(np.median(off)),
        "iid_maximum_check": {"training_realizations": n, "independent_validation_realizations": n,
                              "training_order_1_based": k, "threshold": threshold,
                              "validation_exceedances": validation,
                              "predictive_95pct_count_interval": list(map(int, prediction)),
                              "predictive_mean_count": n * (n + 1 - k) / (n + 1)},
        "stress_test": injection_tests(child_seeds[-1]),
        "global_maxima_per_realization_in_shift_order": all_maxima,
    }
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    main(parser.parse_args().output)
