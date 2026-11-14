"""Replay the two M31-field measurements from the frozen public archive returns."""
import csv
import hashlib
import json
from datetime import date, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
GRID = np.arange(130.0, 1500.01, 0.5)


def read_csv(name, skip=0):
    lines = (DATA / name).read_text().splitlines()[skip:]
    return list(csv.DictReader(lines))


def nightly(rows):
    """One point per integer MJD. Use a floor for calibration and source scatter."""
    data = np.array(rows)
    groups = np.unique(data[:, 0].astype(int))
    data = [np.median(data[data[:, 0].astype(int) == day], axis=0) for day in groups]
    t, magnitude, error = np.array(data).T
    return np.array([t, magnitude, np.hypot(error, 0.10)])


def ztf(year, band):
    return nightly([
        (float(row["mjd"]), float(row["mag"]), float(row["magerr"]))
        for row in read_csv(f"ztf-{year}.csv")
        if row["filtercode"] == band and row["catflags"] == "0" and float(row["magerr"]) < 0.30
    ])


def gaia_dr3():
    return nightly([
        (float(row["g_transit_time"]) + 55197, float(row["g_transit_mag"]),
         1.086 / float(row["g_transit_flux_over_error"]))
        for row in read_csv("gaia-2015.csv")
        if row["g_transit_mag"] and row["variability_flag_g_reject"] == "false"
        and row["rejected_by_photometry"] == "false"
    ])


def gaia_alert():
    rows = []
    for row in read_csv("gaia-2016.csv", skip=1):
        if row["averagemag"] in ("null", "untrusted"):
            continue
        rows.append((float(row["JD(TCB)"]) - 2400000.5, float(row["averagemag"]), 0.20))
    return nightly(rows)


def ps1(band):
    rows = []
    for row in read_csv("ps1-2016.csv"):
        flux, error = float(row["psfFlux"]), float(row["psfFluxErr"])
        if row["filterID"] == band and float(row["psfQfPerfect"]) >= 0.8 and flux > 5 * error:
            rows.append((float(row["obsTime"]), -2.5 * np.log10(flux / 3631), 1.086 * error / flux))
    return nightly(rows)


def wise(year):
    lines = (DATA / f"wise-{year}.tbl").read_text().splitlines()
    header = next(line for line in lines if line.startswith("|"))
    columns = [name.strip() for name in header.strip("|").split("|")]
    rows = []
    for line in lines:
        if not line.strip() or not line.lstrip()[0].isdigit():
            continue
        row = dict(zip(columns, line.split()))
        if (row["cc_flags"][0] == "0" and int(row["qual_frame"]) > 0
                and float(row["saa_sep"]) > 0 and row["moon_masked"][0] == "0"
                and row["w1sigmpro"] != "null" and float(row["w1sigmpro"]) <= 0.25):
            rows.append((float(row["mjd"]), float(row["w1mpro"]), float(row["w1sigmpro"])))
    data = np.array(sorted(rows))
    groups = np.split(data, np.flatnonzero(np.diff(data[:, 0]) > 10) + 1)
    visits = []
    for group in groups:
        if len(group) < 2:
            continue
        t, m, e = group.T
        err = max(0.06, np.median(e) / np.sqrt(len(e)),
                  1.4826 * np.median(np.abs(m - np.median(m))) / np.sqrt(len(e)))
        visits.append((np.median(t), np.median(m), err))
    return np.array(visits).T


def fit(curve, period):
    t, m, e = curve
    angle = (t - 58000) * 2 * np.pi / period
    design = np.column_stack([np.ones(len(t)), np.cos(angle), np.sin(angle)])
    coeff = np.linalg.lstsq(design / e[:, None], m / e, rcond=None)[0]
    residual = m - design @ coeff
    return float(np.sum((residual / e) ** 2)), coeff


def scan(curve):
    return np.array([fit(curve, p)[0] for p in GRID])


def minima(scores, count=3):
    local = np.where((scores[1:-1] < scores[:-2]) & (scores[1:-1] < scores[2:]))[0] + 1
    inds = sorted([0, len(scores) - 1, *local], key=lambda i: scores[i])[:count]
    return [[float(GRID[i]), round(float(scores[i]), 2)] for i in inds]


def phase_max(coeff, period):
    return 58000 + period * ((np.arctan2(coeff[2], coeff[1]) + np.pi) % (2 * np.pi)) / (2 * np.pi)


def main():
    for name, record in json.loads((HERE / "sources.json").read_text()).items():
        if hashlib.sha256((DATA / name).read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"Source digest mismatch for {name}")
    curve = ztf("2015", "zi")
    score = scan(curve)
    period = float(GRID[np.argmin(score)])
    _, coeff = fit(curve, period)
    maximum = phase_max(coeff, period)
    years = np.array([(date(1858, 11, 17) + timedelta(days=float(t))).year for t in curve[0]])
    # Use every represented year with the full period window, including sparsely sampled ones.
    omit = {str(y): minima(scan(curve[:, years != y]), 1)[0][0] for y in np.unique(years)}
    t, m, _ = gaia_dr3()
    infrared = wise("2015")
    sets = {"PS1 z": ps1("4"), "PS1 y": ps1("5"), "Gaia G": gaia_alert(),
            "ZTF i": ztf("2016", "zi"), "ZTF r": ztf("2016", "zr")}
    scores = {key: scan(values) for key, values in sets.items()}
    total = np.sum(list(scores.values()), axis=0)
    result = {
        "2015-07c": {
            "ztf_nights": len(curve[0]), "ztf_best_periods_d_and_scores": minima(score),
            "ztf_peak_to_peak_mag": round(float(2 * np.hypot(coeff[1], coeff[2])), 3),
            "ztf_maximum_mjd": round(float(maximum), 2),
            "discovery_phase_after_max_cycles": round(float(((57216.005-maximum)/period+.5)%1-.5), 3),
            "omit_calendar_year_best_period_d": omit, "wise_visits": len(infrared[0]),
            "wise_best_periods_d_and_scores": minima(scan(infrared)),
            "gaia_dates_mjd_and_g": np.round(np.column_stack([t, m]), 3).tolist(),
        },
        "2016-06a_field": {
            "nights": {k: len(v[0]) for k, v in sets.items()},
            "joint_optical_best_periods_d_and_scores": minima(total),
            "ps1_only_best_periods_d_and_scores": minima(scores["PS1 z"]+scores["PS1 y"]),
            "ztf_gaia_best_periods_d_and_scores": minima(scores["Gaia G"]+scores["ZTF i"]+scores["ZTF r"]),
            "omit_passband_best_period_d": {k: minima(total-v, 1)[0][0] for k, v in scores.items()},
            "wise_visits": len(wise("2016")[0]),
            "wise_best_periods_d_and_scores": minima(scan(wise("2016"))),
        },
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
