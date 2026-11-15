"""Fit individual emission-line cores in the frozen, quality-masked spectral windows."""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


HERE = Path(__file__).resolve().parent
C_KMS = 299792.458
LINES = (4921.931, 5015.678)


def read_pixels():
    rows = np.genfromtxt(HERE / "line-pixels.csv", delimiter=",", names=True,
                         dtype=None, encoding="ascii")
    spectra = {}
    for source in ("XSHOOTER", "BOSS"):
        for line in LINES:
            sel = rows[(rows["source"] == source) & (rows["line_air_A"] == line)]
            epochs = [sel[sel["epoch"] == i] for i in np.unique(sel["epoch"])]
            wavelengths = epochs[0]["wavelength_air_A"]
            if any(not np.array_equal(r["wavelength_air_A"], wavelengths) for r in epochs):
                raise ValueError("Wavelength sampling differs within an instrument.")
            velocity = C_KMS * np.log(wavelengths / line)
            spectra[source, line] = (velocity, np.stack([r["flux"] for r in epochs]),
                                     np.stack([r["ivar"] for r in epochs]))
    return spectra


def fit_core(pixels, epochs, half_width, degree):
    velocity, flux, weight = pixels
    total_weight = weight[epochs].sum(axis=0)
    use = (np.abs(velocity) < half_width) & (total_weight > 0)
    v = velocity[use]
    y = (flux[epochs] * weight[epochs]).sum(axis=0)[use] / total_weight[use]
    precision = np.sqrt(total_weight[use])

    def residual(p):
        amplitude, centre, width = p[:3]
        gaussian = amplitude * np.exp(-0.5 * ((v - centre) / width) ** 2)
        background = np.polynomial.polynomial.polyval(v / 200, p[3:])
        return (gaussian + background - y) * precision

    initial = [max(float(y.max() - np.median(y)), 0.1), 50, 50, np.median(y)] + [0] * degree
    limits = ([0, -130, 10] + [-np.inf] * (degree + 1),
              [np.inf, 190, 175] + [np.inf] * (degree + 1))
    result = least_squares(residual, initial, bounds=limits)
    if not result.success:
        raise RuntimeError(result.message)
    covariance = np.linalg.pinv(result.jac.T @ result.jac)
    return result.x[1], np.sqrt(covariance[1, 1]), result.x[2]


def evaluate(spectra, source, epochs, half_width, degree):
    measurements = [fit_core(spectra[source, line], epochs, half_width, degree)
                    for line in LINES]
    (v1, e1, s1), (v2, e2, s2) = measurements
    return dict(v4922=v1, v5016=v2, sigma4922=s1, sigma5016=s2,
                delta_kms=v1-v2, formal_error_kms=np.hypot(e1, e2))


def main():
    spectra = read_pixels()
    rows = []
    for source, half_width in [("XSHOOTER", 245), ("BOSS", 340)]:
        n = len(spectra[source, LINES[0]][1])
        ix = np.arange(n)
        samples = [("all", ix), ("alternate_A", ix[::2]), ("alternate_B", ix[1::2]),
                   ("first", ix[:n//2]), ("last", ix[n//2:])]
        for label, epochs in samples:
            rows.append(dict(source=source, split=label, count=len(epochs),
                             half_width=half_width, degree=1,
                             **evaluate(spectra, source, epochs, half_width, 1)))
        if source == "XSHOOTER":
            for width in (190, 245, 300, 350):
                for degree in range(3):
                    rows.append(dict(source=source, split="sensitivity", count=n,
                                     half_width=width, degree=degree,
                                     **evaluate(spectra, source, ix, width, degree)))
            individual = [evaluate(spectra, source, [i], half_width, 1)["delta_kms"] for i in ix]
            omitted = [evaluate(spectra, source, np.delete(ix,i), half_width,1)["delta_kms"]
                       for i in ix]
            jackknife_se = np.sqrt((n-1)/n * np.sum((omitted-np.mean(omitted))**2))
            print(json.dumps(dict(positive_individual=int(np.sum(np.array(individual)>0)),
                                  individual_count=n, exposure_jackknife_se_kms=round(jackknife_se, 3))))
    with (HERE / "fits.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows({k: (round(float(v), 4) if isinstance(v, np.floating) else v)
                          for k, v in row.items()} for row in rows)


if __name__ == "__main__":
    main()
