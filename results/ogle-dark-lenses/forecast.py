"""Posterior-predictive I magnitudes at future dates for one event.

Usage: python forecast.py DATADIR EVENT CHAINS.npz OUT.csv

For each saved posterior sample the source and blend fluxes are re-solved on
the photometry, then the model is evaluated at the dates below. The
no-parallax row is the best fit without parallax.
"""
import csv
import pathlib
import sys

import numpy as np

import candidate
import lens

DATES = {"2027-02-15": 11451.5, "2027-04-01": 11496.5, "2027-06-01": 11557.5, "2027-08-01": 11618.5,
         "2027-10-01": 11679.5, "2028-03-01": 11831.5, "2028-06-01": 11923.5, "2028-09-01": 12015.5}


def main():
    data = candidate.prepare(pathlib.Path(sys.argv[1]), sys.argv[2])
    t, flux, sigma, par = data["t"], data["flux"], data["sigma"], data["par"]
    when = np.array(list(DATES.values()))
    future = lens.Parallax(when, data["ra"], data["dec"], data["t0par"])
    chains = np.load(sys.argv[3])
    rows = []
    base = data["base"][0]
    _, fs, fb = lens.model(base, t, flux, sigma)
    tau, beta = lens.trajectory(base, when)
    plain = lens.ZP - 2.5 * np.log10(fs * lens.magnification(tau ** 2 + beta ** 2) + fb)
    baseline = {"no_parallax": lens.ZP - 2.5 * np.log10(fs + fb)}
    predictions = {"no_parallax": np.vstack([plain, plain, plain])}
    for name in ("plus", "minus"):
        sample = chains[name][::4]
        mags = np.empty((len(sample), len(when)))
        level = np.empty(len(sample))
        for i, theta in enumerate(sample):
            _, fs, fb = lens.model(theta, t, flux, sigma, par)
            tau, beta = lens.trajectory(theta, when, future)
            mags[i] = lens.ZP - 2.5 * np.log10(fs * lens.magnification(tau ** 2 + beta ** 2) + fb)
            level[i] = lens.ZP - 2.5 * np.log10(fs + fb)
        predictions[f"u0_{name}"] = np.quantile(mags, [0.5, 0.025, 0.975], axis=0)
        baseline[f"u0_{name}"] = float(np.median(level))
    with open(sys.argv[4], "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["date", "hjd_minus_2450000", "model", "i_median", "i_2.5pct", "i_97.5pct", "baseline_i"])
        for j, (date, hjd) in enumerate(DATES.items()):
            for name, q in predictions.items():
                writer.writerow([date, hjd, name, f"{q[0, j]:.3f}", f"{q[1, j]:.3f}", f"{q[2, j]:.3f}",
                                 f"{baseline[name]:.3f}"])
                rows.append((date, name, q[0, j], q[1, j], q[2, j]))
    for row in rows:
        print("%s %-12s I = %.3f  [%.3f, %.3f]" % row)


if __name__ == "__main__":
    main()
