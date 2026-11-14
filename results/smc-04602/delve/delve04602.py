"""Check OGLE eclipse phases with spatially resolved DELVE-MC photometry."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

PERIOD = 1.7295874
EPOCH = 2460001.4679
WINDOWS = [(2456000, 2457300), (2457300, 2458600), (2458600, 2459700)]
TARGET = "0154m727.1965"
NEIGHBOR = "0154m727.1896"


def load_csv(path):
    with path.open() as source:
        return list(csv.DictReader(source))


def wrap(phase):
    return (phase + 0.5) % 1 - 0.5


def orbital_phase(hjd):
    return ((hjd - EPOCH) / PERIOD) % 1


def mid_hjd(mjd_utc, exposure_seconds):
    """Heliocentric midpoint from a low-precision Earth ephemeris, target ICRS position."""
    jd = mjd_utc + 2400000.5 + exposure_seconds / 172800
    days = jd - 2451545
    longitude = np.deg2rad((280.460 + 0.9856474 * days) % 360)
    anomaly = np.deg2rad((357.528 + 0.9856003 * days) % 360)
    longitude += np.deg2rad(1.915 * np.sin(anomaly) + 0.020 * np.sin(2 * anomaly))
    obliquity = np.deg2rad(23.439 - 0.0000004 * days)
    radius = 1.00014 - 0.01671 * np.cos(anomaly) - 0.00014 * np.cos(2 * anomaly)
    ra, dec = np.deg2rad([15.7324166666667, -72.8732777777778])
    projection = (
        np.cos(longitude) * np.cos(ra) * np.cos(dec)
        + np.sin(longitude) * np.cos(obliquity) * np.sin(ra) * np.cos(dec)
        + np.sin(longitude) * np.sin(obliquity) * np.sin(dec)
    )
    return jd - radius * projection * 499.005 / 86400


def fit_ogle(points):
    time, mag, error = points.T
    phase = orbital_phase(time)
    sigma = np.hypot(error, 0.012)
    baseline = np.percentile(mag, 35)

    def residuals(parameters):
        base, harmonic, dep1, cen1, wid1, dep2, cen2, wid2 = parameters
        prediction = base + harmonic * np.cos(4 * np.pi * phase)
        prediction += dep1 * np.exp(-0.5 * (wrap(phase - cen1) / wid1) ** 2)
        prediction += dep2 * np.exp(-0.5 * (wrap(phase - cen2) / wid2) ** 2)
        return (prediction - mag) / sigma

    bounds = (
        [baseline - 0.5, -0.25, 0.01, -0.22, 0.008, 0.01, 0.28, 0.008],
        [baseline + 0.5, 0.25, 1.8, 0.22, 0.12, 1.8, 0.72, 0.12],
    )
    starts = [
        [baseline, 0.01, 0.45, cen1, 0.03, 0.3, cen2, 0.03]
        for cen1 in (-0.12, 0.0, 0.12) for cen2 in (0.4, 0.53, 0.65)
    ]
    fits = [
        least_squares(residuals, start, bounds=bounds, loss="soft_l1", f_scale=1.8)
        for start in starts
    ]
    return min(fits, key=lambda result: result.cost).x


def unique_exposures(rows, source_id):
    result = {}
    for row in rows:
        if row["objid"] != source_id or row["forced"] != "0" or row["filter"] not in "griz":
            continue
        exposure = row["id"].split("_")[0]
        error = float(row["err"])
        if not 0 < error <= 0.10 or not 0 < float(row["mag"]) < 30:
            continue
        if exposure not in result or error < float(result[exposure]["err"]):
            result[exposure] = row
    return result


def subtract_off_eclipse(rows):
    for band in sorted({row["band"] for row in rows}):
        selected = [row for row in rows if row["band"] == band]
        comparison = [row for row in selected if row["off"]]
        nights = sorted({row["night"] for row in comparison})
        for field in ("m", "ctrl", "delta"):
            nightly = [
                np.median([row[field] for row in comparison if row["night"] == night])
                for night in nights
            ]
            baseline = np.median(nightly) if nightly else np.nan
            for row in selected:
                row[field + "_excess"] = float(row[field] - baseline)


def analyze(directory, ogle):
    observations = np.loadtxt(ogle)
    observations = observations[(observations[:, 2] > 0) & (observations[:, 2] <= 0.10)]
    fits = [
        fit_ogle(observations[(observations[:, 0] >= start) & (observations[:, 0] < end)])
        for start, end in WINDOWS
    ]
    target = unique_exposures(load_csv(directory / "delve-meas.csv"), TARGET)
    control = unique_exposures(load_csv(directory / "delve-control.csv"), NEIGHBOR)
    exptime = {
        record["expnum"]: float(record["exptime"])
        for file in sorted(directory.glob("delve-expo*.csv")) for record in load_csv(file)
    }
    rows = []
    for exposure in sorted(target.keys() & control.keys() & exptime.keys()):
        star, neighbor = target[exposure], control[exposure]
        time = mid_hjd(float(star["mjd"]), exptime[exposure])
        phase = orbital_phase(time)
        window = next(i for i, (start, end) in enumerate(WINDOWS) if start <= time < end)
        distances = np.abs(wrap(phase - fits[window][[3, 6]]))
        near = int(np.argmin(distances))
        mag, ctrl = float(star["mag"]), float(neighbor["mag"])
        rows.append(dict(
            exp=exposure, hjd=time, window=window, band=star["filter"], phase=phase,
            core=near + 1 if distances[near] < 0.03 else 0, off=bool(min(distances) > 0.12),
            dphase=float(distances[near]), m=mag, ctrl=ctrl, err=float(star["err"]),
            ctrl_err=float(neighbor["err"]), night=int(time), delta=mag - ctrl,
        ))
    subtract_off_eclipse(rows)
    result = {
        "OGLE_windows": WINDOWS,
        "OGLE_parameters": np.array(fits).tolist(),
        "unique_matched_exposures": len(rows),
        "distinct_nights": len({row["night"] for row in rows}),
        "core": {},
    }
    for eclipse in (1, 2):
        core = [row for row in rows if row["core"] == eclipse]
        records = []
        for night in sorted({row["night"] for row in core}):
            sample = [row for row in core if row["night"] == night]
            excess = {
                field: round(float(np.median([row[field + "_excess"] for row in sample])), 4)
                for field in ("m", "ctrl", "delta")
            }
            records.append(dict(
                hjd_night=night, n=len(sample), bands="".join(sorted({r["band"] for r in sample})),
                excess=excess, min_distance=min(r["dphase"] for r in sample),
            ))
        result["core"][eclipse] = records
    return result, rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("ogle_iv_I", type=Path)
    args = parser.parse_args()
    result, rows = analyze(args.directory, args.ogle_iv_I)
    (args.directory / "delve-analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    with (args.directory / "delve-analyzed.csv").open("w") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} matched exposures on {result['distinct_nights']} nights")
