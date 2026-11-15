"""Compare a common intrinsic gap with a simple frequency-dependent scattering tail."""
import argparse
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import exponnorm

from replay import HERE, inputs, save_csv


def shape(t, width, tau):
    if tau < 0.002:
        return np.exp(-0.5 * (t / width) ** 2)
    convolved = exponnorm.pdf(t, tau / width, loc=0, scale=width)
    return convolved / convolved.max()


def spectra(cache):
    data = []
    for label in [85, 119, 154, 185, 215]:
        flux = np.loadtxt(cache / f"{label}MHz_dedisper.txt")
        nu = np.loadtxt(cache / f"{label}MHz_freq.txt") * 1e-6
        time = np.arange(flux.shape[1]) * 0.5
        on = (time >= 32) & (time <= 55)
        off = (time >= 56) & (time <= 63)
        assert np.all(flux[:, on | off] != 0)
        avg = flux.mean(axis=0)
        data.append((nu.mean(), time[on], avg[on], np.std(avg[off])))
    return data


def compare(data, model):
    n = len(data)
    # Each band gets its own leading location, two heights, two widths,
    # baseline and slope. The frequency law applies only to scattering.
    nuisance = [40, 12, 9, 2.3, 2.3, 0, 0]
    lower = list(np.tile([34, .01, .01, .2, .2, -15, -4], n))
    upper = list(np.tile([43, 80, 80, 8, 8, 15, 4], n))
    if model == "separate_gaps":
        lower += [3] * n
        upper += [12] * n
        starts = [[7] * n]
    elif model == "one_gap":
        lower += [3]
        upper += [12]
        starts = [[7]]
    else:
        lower += [3, 0]
        upper += [12, 15]
        starts = [[7, tau] for tau in [0, 0.1, 0.5, 1.5, 3, 6, 12]]

    def residual(p):
        pieces = []
        for i,(nu,t,y,err) in enumerate(data):
            t0, a, b, w1, w2, baseline, slope = p[7*i:7*i+7]
            gap = p[7*n+i] if model == "separate_gaps" else p[7*n]
            tau = p[-1]*(nu/88)**-4 if model == "one_gap_and_tail" else 0
            fit = a*shape(t-t0,w1,tau) + b*shape(t-t0-gap,w2,tau) + baseline + slope*(t-43)
            pieces.append((fit - y) / err)
        return np.concatenate(pieces)

    fits = [least_squares(residual, nuisance*n + start, bounds=(lower,upper), max_nfev=350)
            for start in starts]
    best = min(fits, key=lambda fit: fit.cost)
    if not best.success:
        raise RuntimeError(best.message)
    return dict(bands=n,model=model,parameters=len(best.x),samples=len(best.fun),
                weighted_rss=best.fun@best.fun,
                gaps_s=" ".join(f"{x:.3f}" for x in best.x[7*n:(8*n if model=="separate_gaps" else 7*n+1)]),
                tau_at_88_s=best.x[-1] if model=="one_gap_and_tail" else 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=HERE)
    args = parser.parse_args()
    inputs(args.cache)
    all_spectra = spectra(args.cache)
    rows = [compare(data, model) for data in [all_spectra, all_spectra[1:]]
            for model in ["one_gap", "one_gap_and_tail", "separate_gaps"]]
    args.out.mkdir(parents=True, exist_ok=True)
    save_csv(args.out / "scattering_models.csv", rows)
