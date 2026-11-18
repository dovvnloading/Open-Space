"""Check spectral holdout predictions conditional on the total counts in each bin."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, xlogy

from analyze import HERE, BASELINE, SPECTRAL, NEGATIVE_SPECTRAL, window, analysis, NULL


def check(spectra, band, sigma, t, interval):
    ref, obs = window(t, BASELINE), window(t, interval)
    reference = spectra[ref].sum(axis=0)
    target = spectra[obs].sum(axis=0)
    bins = np.flatnonzero((band >= 1200) & (band <= 1800))
    rows = []
    for group in range(len(bins) // 16):
        selected = bins[group * 16:group * 16 + 16]
        r = reference[selected]
        if r.sum() >= 40:
            rows.append((group, float(np.average(band[selected], weights=r)), int(r.sum()),
                         int(target[selected].sum()), float(np.average(sigma[selected], weights=r) * 1e16)))
    g, wave, low, high, absorb = np.asarray(rows).T
    total = low + high
    train = (g % 2) == 0
    scale = np.log(obs.sum() / ref.sum())
    a = np.column_stack((np.ones(len(g)), (wave - 1500) / 300))
    fitted = {}
    for name, design, bounds in (
        ("gray_color", a, [(-5, 5), (-3, 3)]),
        ("gray_color_water", np.column_stack((a, -absorb)), [(-5, 5), (-3, 3), (0, 20)]),
    ):
        def objective(params):
            eta = scale + design[train] @ params
            return (float(np.sum(total[train] * np.logaddexp(0, eta) - high[train] * eta)),
                    design[train].T @ (total[train] * expit(eta) - high[train]))

        fit = minimize(objective, np.zeros(design.shape[1]), method="L-BFGS-B", jac=True, bounds=bounds)
        p = expit(scale + design[~train] @ fit.x)
        y, n = high[~train], total[~train]
        deviance = 2 * np.sum(xlogy(y, y / (n * p)) + xlogy(n-y, (n-y) / (n * (1-p))))
        fitted[name] = {"params": list(map(float, fit.x)), "converged": bool(fit.success),
                        "heldout_deviance": float(deviance)}
    return {"training_blocks": int(train.sum()), "heldout_blocks": int((~train).sum()), "models": fitted,
            "heldout_improvement_with_water": fitted["gray_color"]["heldout_deviance"]
            - fitted["gray_color_water"]["heldout_deviance"]}


def replay(root):
    result, time_series, cube = analysis(root)
    core = np.array(result["array_support"]["core_spatial_rows_zero_based"]) - 30
    spectra = np.where(cube != NULL, cube, 0).astype(np.int64)
    band = np.loadtxt(HERE / "calibration/fuv-wavelength-bin-centers.csv", delimiter=",").ravel()
    water = np.loadtxt(HERE / "calibration/H2O-Mota-2005-298K.txt")
    sigma = np.interp(band / 10, water[:, 0], water[:, 1], left=np.nan, right=np.nan)
    t = time_series[:, 0]
    return {"encounter": check(spectra.sum(axis=1), band, sigma, t, SPECTRAL),
            "encounter_core_rows": check(spectra[:, core, :].sum(axis=1), band, sigma, t, SPECTRAL),
            "quiet": check(spectra.sum(axis=1), band, sigma, t, NEGATIVE_SPECTRAL)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = replay(args.raw)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print({key: entry["heldout_improvement_with_water"] for key, entry in result.items()})


if __name__ == "__main__":
    main()
