"""Replay KELT-8 transit times using the resolved Gaia neighbour as a template."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
from functools import partial
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.time import Time
from scipy.optimize import minimize_scalar

from pixels import locate


HERE = Path(__file__).resolve().parent
REFERENCE_BTJD = 986.46737  # Kokori et al. (2023), BJD_TDB minus 2457000
PLANET_DAYS = 3.24408156
DEFAULT_SHAPE = (0.1444, 0.0277)  # total and ingress duration, days (Fulton et al.)
EPOCH = "orbit"
EVENT_FIELDS = ["sector", EPOCH, "raw_min", "corrected_min", "predicted_bias_min",
                "raw_formal_err_min", "corrected_formal_err_min", "raw_depth_ppm",
                "corrected_depth_ppm"]


def read_inputs(cache, fetch):
    paths = {}
    for entry in json.loads((HERE / "inputs.json").read_text()):
        name = entry["name"]
        path = (cache if name.endswith(".fits") else HERE) / name
        if not path.exists() and fetch and name.endswith(".fits"):
            with urlopen(entry["url"], timeout=120) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != entry["sha256"]:
                raise ValueError(f"Remote product has changed: {name}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"Input hash mismatch: {path}")
        paths[name] = path
    return paths


def relative_binary_flux(btjd, model, phase_shift=0.0):
    """Evaluate Gaia's two-Gaussian + ellipsoidal G-band model at TESS times."""
    reference_tcb = 2455197.5 + float(model["reference_time"])
    tcb = Time(2457000.0 + btjd, format="jd", scale="tdb").tcb.jd
    phase = (tcb - reference_tcb) * float(model["frequency"]) + phase_shift
    mag = np.zeros_like(btjd)
    for i in (1, 2):
        prefix = f"geom_model_gaussian{i}_"
        distance = (phase - float(model[prefix + "phase"]) + 0.5) % 1.0 - 0.5
        mag += float(model[prefix + "depth"]) * np.exp(
            -0.5 * (distance / float(model[prefix + "sigma"])) ** 2
        )
    ellipsoidal = "geom_model_cosine_half_period_"
    mag += float(model[ellipsoidal + "amplitude"]) * np.cos(
        4 * np.pi * (phase - float(model[ellipsoidal + "phase"]))
    )
    return 10 ** (-0.4 * mag)


def sector_dilution(time, flux, error, binary, outside):
    """Regress out-of-transit data against the Gaia template after local detrending."""
    target, template, errors = [], [], []
    day = np.floor(time)
    for date in np.unique(day):
        use = outside & (day == date)
        if use.sum() < 100:
            continue
        design = np.column_stack((np.ones(use.sum()), time[use] - date))
        for series, output in [(flux, target), (binary, template)]:
            coef = np.linalg.lstsq(design / error[use, None], series[use] / error[use],
                                   rcond=None)[0]
            output.extend(series[use] - design @ coef)
        errors.extend(error[use])
    y, x, sigma = map(np.asarray, [target, template, errors])
    keep = np.ones(len(y), dtype=bool)
    for _ in range(4):
        scale = np.sum(y[keep] * x[keep] / sigma[keep]**2) / np.sum((x[keep] / sigma[keep])**2)
        resid = y - scale * x
        dist = abs(resid - np.median(resid))
        keep = dist < 5 * 1.4826 * np.median(dist)
    return scale, np.corrcoef(y, x)[0, 1], len(y)


def trapezoid(offset, shape):
    duration, ingress = shape
    return np.clip((duration / 2 - np.abs(offset)) / ingress, 0.0, 1.0)


def fit_center(offset, flux, error, shape):
    def profile(shift):
        design = np.column_stack((np.ones(len(offset)), offset, trapezoid(offset - shift, shape)))
        coef = np.linalg.lstsq(design / error[:, None], flux / error, rcond=None)[0]
        return (flux - design @ coef) / error, -coef[-1]

    opt = minimize_scalar(lambda shift: np.sum(profile(shift)[0]**2),
                          bounds=(-0.018, 0.018), method="bounded", options={"xatol": 1e-8})
    shift = opt.x
    step = 1e-5
    derivative = (profile(shift + step)[0] - profile(shift - step)[0]) / (2 * step)
    return shift * 1440, 1440 / np.linalg.norm(derivative), profile(shift)[1] * 1e6


def measure_sector(path, model, column="PDCSAP_FLUX", phase_shift=0., span=0.23,
                   shape=DEFAULT_SHAPE):
    with fits.open(path) as product:
        sector = product[0].header["SECTOR"]
        data = product[1].data
        use = ((data["QUALITY"] == 0) & np.isfinite(data["TIME"]) & np.isfinite(data[column])
               & np.isfinite(data[column + "_ERR"]) & (data[column + "_ERR"] > 0))
        time, flux, error = (data[k][use].astype(float)
                             for k in ("TIME", column, column + "_ERR"))
    median = np.median(flux)
    flux, error = flux / median, error / median
    orbit = np.rint((time - REFERENCE_BTJD) / PLANET_DAYS).astype(int)
    offset = time - REFERENCE_BTJD - orbit * PLANET_DAYS
    binary = relative_binary_flux(time, model, phase_shift)
    scale, corr, n_calibration = sector_dilution(time, flux, error, binary, abs(offset) > 0.13)
    rows = []
    cut = np.array([-span, -0.085, -0.078, -0.04, 0.04, 0.078, 0.085, span])
    regions = [(cut[i], cut[j]) for i, j in [(0, 1), (2, 3), (3, 4), (4, 5), (6, 7)]]
    for n in np.unique(orbit):
        selected = (orbit == n) & (abs(offset) < span)
        t, y, err = offset[selected], flux[selected], error[selected]
        if any(np.sum((t > lo) & (t < hi)) < 10 for lo, hi in regions):
            continue
        additive = scale * binary[selected]
        raw = fit_center(t, y, err, shape)
        corrected = fit_center(t, y - additive, err, shape)
        injection = 1 - 0.013 * trapezoid(t, shape) + additive
        bias = fit_center(t, injection, err, shape)[0]
        values = (sector, int(n), raw[0], corrected[0], bias, raw[1], corrected[1],
                  raw[2], corrected[2])
        rows.append(dict(zip(EVENT_FIELDS, values)))
    return rows, dict(sector=sector, cadences=len(time), calibration_rows=n_calibration,
                      fitted_scale=scale, out_of_transit_r=corr, transits=len(rows))


def detrend(epochs, values):
    design = np.column_stack((np.ones(len(epochs)), epochs - np.mean(epochs)))
    return values - design @ np.linalg.lstsq(design, values, rcond=None)[0]


def evaluate(rows):
    summary = {}
    for group, keep in [("catalog_sectors", lambda s: s <= 54),
                        ("later_sectors", lambda s: s >= 80),
                        ("all", lambda s: True)]:
        sample = [row for row in rows if keep(row["sector"])]
        epoch = np.array([row[EPOCH] for row in sample])
        raw, clean, pred = (detrend(epoch, np.array([row[k] for row in sample]))
                            for k in ["raw_min", "corrected_min", "predicted_bias_min"])
        summary[group] = dict(n=len(sample), raw_rms_min=np.std(raw),
                              corrected_rms_min=np.std(clean), correlation=np.corrcoef(raw, pred)[0, 1],
                              slope=raw @ pred / (pred @ pred), phase_agreement=int(np.sum(raw * pred > 0)),
                              prediction_residual_rms_min=np.std(raw - pred))
    return summary


def catalog_matches(paths, rows):
    old = {r[EPOCH]: r for r in rows if r["sector"] <= 54}
    files = {"wang_times": ("wang-times.csv", ",", 1), "wang_espinoza": ("wang-espinoza.dat", None, 0)}
    stats = {}
    for key, (name, delimiter, index) in files.items():
        time = np.loadtxt(paths[name], delimiter=delimiter)[:, index]
        epoch = np.rint((time - 2457000 - REFERENCE_BTJD) / PLANET_DAYS).astype(int)
        pred = detrend(epoch, np.array([old[n]["predicted_bias_min"] for n in epoch]))
        obs = detrend(epoch, (time - 2457000 - REFERENCE_BTJD - epoch * PLANET_DAYS) * 1440)
        slope = pred @ obs / (pred @ pred)
        stats[key] = dict(n=len(epoch), correlation=np.corrcoef(pred, obs)[0, 1], slope=slope,
                          rms_min=np.std(obs), residual_rms_min=np.std(obs - slope * pred))
    return stats


def write_csv(path, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path, help="directory for frozen NASA FITS")
    parser.add_argument("--fetch", action="store_true", help="retrieve missing FITS after checking hashes")
    parser.add_argument("--out", type=Path, required=True, help="output directory")
    args = parser.parse_args()
    paths = read_inputs(args.cache, args.fetch)
    with paths["gaia-eb.csv"].open() as f:
        model = next(csv.DictReader(f))
    products = [path for name, path in paths.items() if name.endswith("_lc.fits")]
    options = {"primary": {}, "sap_flux": {"column": "SAP_FLUX"}, "wider_window": {"span": 0.30},
               "longer_ingress": {"shape": (0.1444, 0.0332)},
               "wrong_eb_phase": {"phase_shift": 0.125}}
    args.out.mkdir(parents=True, exist_ok=True)
    summary = {}
    for label, config in options.items():
        rows, sectors = [], []
        for path in products:
            measured, info = measure_sector(path, model, **config)
            rows.extend(measured)
            sectors.append(info)
        summary[label] = evaluate(rows)
        if label == "primary":
            write_csv(args.out / "events.csv", rows)
            write_csv(args.out / "sectors.csv", sectors)
            summary["published_times"] = catalog_matches(paths, rows)
    positions = []
    for name, path in paths.items():
        if name.endswith("_tp.fits"):
            positions.extend(locate(path, paths["field.csv"], partial(relative_binary_flux, model=model),
                                    REFERENCE_BTJD, PLANET_DAYS, partial(trapezoid, shape=DEFAULT_SHAPE)))
    write_csv(args.out / "pixel-centroids.csv", positions)
    rounded = {k: {s: {stat: round(float(v), 5) if not isinstance(v, int) else v
                       for stat, v in stats.items()} for s, stats in groups.items()}
               for k, groups in summary.items()}
    (args.out / "summary.json").write_text(json.dumps(rounded, indent=2) + "\n")


if __name__ == "__main__":
    main()
