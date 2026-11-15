"""Fixed-redshift aperture equivalent widths in independent archival spectra."""
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
GROUPS = {
    "MgII": [2796.3543, 2803.5315],
    "FeIIa": [2344.2139],
    "FeIIb": [2374.4603, 2382.7642],
    "FeIIc": [2586.6494, 2600.1725],
    "MgI": [2852.9631],
}


def measure(wave, flux, ivar, mask, z, group, half_width=210, continuum="linear"):
    centers = np.asarray(GROUPS[group]) * (1 + z)
    x = (wave - centers.mean()) / 40
    design = np.column_stack([np.ones_like(x), x]) if continuum == "linear" else np.ones((len(x), 1))
    # Omit all measured transitions, even when a second fitting group overlaps this continuum.
    all_centers = np.array([line for lines in GROUPS.values() for line in lines]) * (1 + z)
    side = (np.abs(x) < 1.0) & (abs(wave[:, None] - all_centers).min(axis=1) > 8)
    good = (ivar > 0) & (mask == 0) & np.isfinite(flux)
    side &= good
    for _ in range(2):
        if side.sum() < 10:
            return []
        normal = design[side].T * ivar[side]
        cov = np.linalg.inv(normal @ design[side])
        model = design @ (cov @ (normal @ flux[side]))
        side &= np.abs(flux - model) * np.sqrt(np.maximum(ivar, 0)) < 4

    edges = np.r_[wave[0] - (wave[1]-wave[0])/2, (wave[1:] + wave[:-1]) / 2,
                  wave[-1] + (wave[-1]-wave[-2])/2]
    output = []
    for rest, center in zip(GROUPS[group], centers):
        half = half_width / 299792.458 * center
        widths = np.maximum(0, np.minimum(edges[1:], center + half) - np.maximum(edges[:-1], center - half))
        line = widths > 0
        if not good[line].all() or (model[line] <= 0).any() or line.sum() < 3:
            continue
        weight = widths[line] / (1+z)
        ew = np.sum(weight * (1 - flux[line] / model[line]))
        grad = (weight * flux[line] / model[line]**2) @ design[line]
        error = np.sqrt(np.sum((weight / model[line])**2 / ivar[line]) + grad @ cov @ grad)
        output.append((rest, ew, error))
    return output


def main():
    visits = {row["spectrum"]: row for row in csv.DictReader((HERE / "visits.csv").open())}
    with (HERE / "line_measurements.csv").open("w", newline="") as stream:
        out = csv.writer(stream, lineterminator="\n")
        out.writerow(["beam","spectrum","line","rest_A","halfwidth_kms","continuum","EW_rest_A","sigma_A","EW_over_sigma"])
        for spec, visit in visits.items():
            s = np.genfromtxt(HERE / "pixels" / f"{spec}.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
            for group in GROUPS:
                for width, model in [(210,"linear"), (150,"linear"), (275,"linear"), (210,"constant")]:
                    rows = measure(s["wave_A"], s["flux"], s["ivar"], s["mask"], float(visit["z_abs_desi"]), group, width, model)
                    for rest, ew, error in rows:
                        out.writerow([visit["beam"], spec, group, rest, width, model, f"{ew:.4f}", f"{error:.4f}", f"{ew/error:.3f}"])


if __name__ == "__main__":
    main()
