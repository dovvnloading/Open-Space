"""Digitize the iMet public panel markers and test a fixed period on held-out pages."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from PIL import Image
from scipy.cluster.hierarchy import fclusterdata
from scipy.ndimage import distance_transform_edt, maximum_filter
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
PAGES = ["1st", "2nd", "3rd", "4th", "5th", "6th"]
COLORS = {"blue": (0, 0, 255), "black": (0, 0, 0), "yellow": (255, 215, 0), "green": (0, 128, 0)}
CHANNEL_RGB = {**COLORS, "red": (255, 0, 0)}
# (flux1, row1, flux2, row2) read at labeled tick centers on 1500 × 844 reference images.
TICKS = [(1, 689, 8, 190), (1, 695, 8, 190), (1, 687, 8, 190),
         (0, 714, 8, 186), (0, 713, 8, 184), (0, 713, 8, 184)]


def digitize(path, panel, color, ticks=TICKS):
    pixels = np.asarray(Image.open(path).convert("RGB"), dtype="int16")
    height, width = pixels.shape[:2]
    left, right = 0.10 * width, 0.85 * width
    mask = np.max(np.abs(pixels - CHANNEL_RGB[color]), axis=2) < 85
    mask[:round(0.15 * height)] = False
    mask[round(0.865 * height):] = False
    mask[:, np.arange(width) < left + 3] = False
    mask[:, np.arange(width) > right - 3] = False
    # Markers have thicker interiors than the connecting segments. Several core
    # maxima on one plotting symbol can arise from JPEG compression.
    depth = distance_transform_edt(mask)
    core = (depth == maximum_filter(depth, size=(9, 5))) & (depth >= 2.8 * width / 1500)
    row, col = np.where(core)
    labels = fclusterdata(np.column_stack((col, row)), 3.8 * width / 1500, criterion="distance")
    centers = np.array([(np.median(col[labels == j]), np.median(row[labels == j]))
                        for j in np.unique(labels)])
    flux1, row1, flux2, row2 = ticks[panel]
    mjd = 56000 + 1000 * panel + (centers[:, 0] - left) * 1000 / (right - left)
    flux = flux1 + (centers[:, 1] * 844 / height - row1) * (flux2 - flux1) / (row2 - row1)
    points = np.column_stack((mjd, flux))
    bins = np.floor(mjd / 5).astype(int)
    return np.array([np.median(points[bins == j], axis=0) for j in np.unique(bins)])


def design(time, period):
    angle = 2 * np.pi * (time - 59000) / period
    return np.column_stack((np.ones_like(angle), np.sin(angle), np.cos(angle),
                            np.sin(2 * angle), np.cos(2 * angle)))


def fit(data, low=90, high=130):
    time, flux = data.T

    def residual(period):
        matrix = design(time, period)
        coeff = np.linalg.lstsq(matrix, flux, rcond=None)[0]
        return np.sum((flux - matrix @ coeff) ** 2)

    grid = np.arange(low, high, 0.05)
    best = min(grid, key=residual)
    period = minimize_scalar(residual, bounds=(best - 0.06, best + 0.06), method="bounded").x
    return period, np.linalg.lstsq(design(time, period), flux, rcond=None)[0]


def serialize(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((ROOT / "panels.json").read_text())
    paths = []
    for page in PAGES:
        filename = "28.02-04_" + page + ".jpg"
        path = args.cache / filename
        if not path.exists():
            with urlopen(manifest[filename]["url"], timeout=45) as response:
                content = response.read()
            if hashlib.sha256(content).hexdigest() != manifest[filename]["sha256"]:
                raise ValueError(f"Live image changed: {filename}. Use the archived measurements.")
            path.write_bytes(content)
        if hashlib.sha256(path.read_bytes()).hexdigest() != manifest[filename]["sha256"]:
            raise ValueError(f"Image digest mismatch: {filename}")
        paths.append(path)

    curves = {color: [digitize(path, i, color) for i, path in enumerate(paths)] for color in COLORS}
    with (args.out / "digitized.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["color", "panel", "mjd_approx", "flux_jy_approx"])
        for color, panels in curves.items():
            for page, data in zip(PAGES, panels):
                writer.writerows((color, page, f"{t:.3f}", f"{v:.4f}") for t, v in data)

    results = {}
    for color, panels in curves.items():
        training = np.vstack(panels[2:5])
        period, coeff = fit(training)
        time = np.linspace(59000, 59000 + period, 5000)
        prediction = design(time, period) @ coeff
        panels_summary = []
        for page, data in zip(PAGES, panels):
            t, y = data.T
            panels_summary.append(dict(page=page, bins=len(t), first_mjd=float(min(t)), last_mjd=float(max(t)),
                                       r_with_frozen_model=float(np.corrcoef(y, design(t, period) @ coeff)[0, 1])))
        results[color] = dict(period_days=period, maximum_mjd=time[np.argmax(prediction)],
                              peak_to_trough_jy=np.ptp(prediction), coefficients=coeff.tolist(), panels=panels_summary)
    # The neighboring velocity traces share the receiver, observation times and calibration.
    # Test each against the blue model, including after differencing log flux from green.
    period = results["blue"]["period_days"]
    coeff = results["blue"]["coefficients"]
    from scipy.spatial import cKDTree
    ratios = []
    for blue, green in zip(curves["blue"], curves["green"]):
        delta, near = cKDTree(green[:, :1]).query(blue[:, :1])
        matched = (delta < 2.5) & (green[near, 1] > 0.25)
        t, y = blue[matched].T
        z = np.log(y / green[near[matched], 1])
        expectation = design(t, period) @ coeff
        ratios.append(dict(bins=len(t), blue_green_log_ratio_r=float(np.corrcoef(z, expectation)[0, 1]),
                           green_r=float(np.corrcoef(green[near[matched], 1], expectation)[0, 1])))
    results["calibration_control"] = ratios
    results["blue_leave_one_training_panel_out_days"] = [
        fit(np.vstack([curves["blue"][j] for j in (2, 3, 4) if j != omit]))[0] for omit in (2, 3, 4)
    ]
    serialize(args.out / "measurements.json", results)


if __name__ == "__main__":
    main()
