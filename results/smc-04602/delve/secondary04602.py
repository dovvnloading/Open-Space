"""Fit centers of the resolved DECam secondary eclipse across passbands and shapes."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from delve04602 import EPOCH, PERIOD


def fit(rows, bands, measurement, shape):
    selected = [
        row for row in rows if row["band"] in bands and int(row["night"]) == 2459151
    ]
    channels = np.array([bands.index(row["band"]) for row in selected])
    phase = np.array([float(row["phase"]) for row in selected])
    mag = np.array([float(row[measurement + "_excess"]) for row in selected])
    variance = np.array([float(row["err"]) ** 2 for row in selected]) + 0.012 ** 2
    if measurement == "delta":
        variance += np.array([float(row["ctrl_err"]) ** 2 for row in selected])
    sigma = np.sqrt(variance)
    nband = len(bands)

    def residuals(parameters):
        center, width = parameters[:2]
        amplitude = parameters[2:2 + nband]
        baseline = parameters[2 + nband:]
        distance = np.abs((phase - center) / width)
        if shape == "triangle":
            profile = np.maximum(0, 1 - distance / 2.5)
        else:
            power = 2 if shape == "gaussian" else 3
            profile = np.exp(-0.5 * distance ** power)
        model = baseline[channels] + amplitude[channels] * profile
        return (mag - model) / sigma

    result = least_squares(
        residuals, [0.61, 0.022] + [0.35] * nband + [0.0] * nband,
        bounds=(
            [0.574, 0.005] + [0.03] * nband + [-0.2] * nband,
            [0.675, 0.09] + [1.0] * nband + [0.2] * nband,
        ), loss="soft_l1",
    )
    center = result.x[0]
    cycle = np.rint((2459151.5 - EPOCH) / PERIOD - center)
    return dict(
        center_phase=center, hjd_utc=EPOCH + (cycle + center) * PERIOD,
        width_scale_phase=result.x[1], n=len(phase), depth_mag=result.x[2:2 + nband].tolist(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    directory = parser.parse_args().directory
    with (directory / "delve-analyzed.csv").open() as file:
        rows = list(csv.DictReader(file))
    fits = {
        measurement: {
            shape: {
                band: fit(rows, band, measurement, shape) for band in ("griz", "g", "r", "i", "z")
            }
            for shape in ("gaussian", "triangle", "supergaussian")
        }
        for measurement in ("m", "delta")
    }
    (directory / "delve-secondary.json").write_text(json.dumps(fits, indent=2) + "\n")
    print(f"Secondary center {fits['m']['gaussian']['griz']['center_phase']:.6f} cycle")
