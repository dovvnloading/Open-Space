#!/usr/bin/env python3
"""Independent NASA OMNI indices vs split-group MeerKAT Faraday residuals.

Reads NASA SPDF OMNI2 hourly ASCII *for the actual observation years*.
OMNI Kp comes from GFZ; F10.7 from NRC; Dst from Kyoto WDC. These are
independent environmental indices, NOT independent ionospheric TEC/RM maps.
Downloading is opt-in; --directory can contain previously fetched files.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / "meerkat-shared-faraday-rotation"
BASE = "https://cdaweb.gsfc.nasa.gov/pub/data/omni/low_res_omni/"

# OMNI2 documentation: 39=Kp*10, 41=Dst nT, 51=F10.7 sfu (1-based).
# https://omniweb.gsfc.nasa.gov/html/ow_data.html
# https://github.com/st-bender/pyspaceweather/blob/master/omnie_format.txt
KPIDX, DSTIDX, FIDX = 38, 40, 50


def fetch_year(year, directory, allow_download=False):
    path = directory / f"omni2_{year}.dat"
    if not path.exists():
        if not allow_download:
            raise FileNotFoundError(f"{path} absent; supply --fetch or --directory")
        directory.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(
            BASE + path.name, headers={"User-Agent": "MeerKAT-common-mode-replay/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response:
            data = response.read()
        if len(data) < 1_000_000:
            raise ValueError(f"Unexpectedly short OMNI file: {path.name} ({len(data)})")
        path.write_bytes(data)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return path, digest


def indices_from_files(paths):
    all_hours = defaultdict(lambda: {"kp": [], "dst": [], "f107": []})
    for path in paths:
        for line in path.read_text(errors="replace").splitlines():
            p = line.split()
            if len(p) < 51:
                continue
            try:
                yr, doy, hour = map(int, p[:3])
                date = dt.date(yr, 1, 1) + dt.timedelta(days=doy - 1)
                kp = float(p[KPIDX]) / 10.0
                dst = float(p[DSTIDX])
                f107 = float(p[FIDX])
            except (ValueError, OverflowError):
                continue
            if not (1 <= doy <= 366 and 0 <= hour < 24):
                continue
            day = all_hours[date]
            if 0 <= kp <= 9:
                day["kp"].append(kp)
            if -1000 < dst < 500:
                day["dst"].append(dst)
            if 40 < f107 < 900:
                day["f107"].append(f107)
    daily = {}
    for date, vals in all_hours.items():
        if all(len(a) >= 1 for a in vals.values()):
            daily[date] = {"kp_max": max(vals["kp"]),
                           "dst_min": min(vals["dst"]),
                           "f107_med": float(np.median(vals["f107"]))}
    return daily


def pearson(x, y):
    return float(np.corrcoef(x, y)[0, 1]) if len(x) >= 3 else float("nan")


def test_external_indices(daily):
    with (PRIOR / "dates.csv").open(newline="") as f:
        dates = list(csv.DictReader(f))
    matched = []
    for r in dates:
        date = dt.date(1858,11,17) + dt.timedelta(days=int(float(r["mjd_utc"])))
        if date not in daily:
            continue
        a, b = float(r["resid_A"]), float(r["resid_B"])
        matched.append({"date": date.isoformat(), "year": date.year,
                        "mjd": int(float(r["mjd_utc"])), "A": a, "B": b,
                        "common": .5*(a+b), "abs_common": abs(.5*(a+b)),
                        **daily[date]})
    if len(matched) != len(dates):
        raise ValueError(f"Missing OMNI indices for {len(dates)-len(matched)} of {len(dates)} days")
    # Prespecified contrasts; no hunting for an optimal lag or weather threshold.
    indices = ("kp_max", "dst_min", "f107_med")
    r = {key: {"signed_r": pearson([z["common"] for z in matched],
                                  [z[key] for z in matched]),
               "absolute_r": pearson([z["abs_common"] for z in matched],
                                    [z[key] for z in matched])}
         for key in indices}
    storm = [z for z in matched if z["kp_max"] >= 5 or z["dst_min"] <= -50]
    quiet = [z for z in matched if z["kp_max"] < 3 and z["dst_min"] > -30]
    # Year-holdout: train 2020-21, evaluate 2022-23. No heldout refit.
    train = [z for z in matched if z["year"] <= 2021]
    test = [z for z in matched if z["year"] >= 2022]
    def fitted(terms, ridge=0):
        X = np.array([[z[k] for k in terms] for z in train])
        mean, std = X.mean(axis=0), X.std(axis=0)
        std = np.where(std > 0, std, 1.)
        Z = np.column_stack([np.ones(len(X)), (X-mean)/std])
        y = np.array([z["B"] for z in train])
        reg = np.eye(len(terms)+1) * ridge
        reg[0,0] = 0.
        b = np.linalg.solve(Z.T@Z+reg, Z.T@y)
        def predict(group):
            xx = np.array([[z[k] for k in terms] for z in group])
            return b[0]+((xx-mean)/std)@b[1:]
        return {"intercept": float(b[0]),
                "betas_z": {k:float(b[i+1]) for i,k in enumerate(terms)}}, predict
    models = {
        "no_correction": lambda g: np.zeros(len(g)),
        "unscaled_peers": lambda g: np.array([z["A"] for z in g])}
    coefficients={}
    for name, terms in [("indices_only",indices),
                        ("peers_only",("A",)),
                        ("peers_plus_indices",("A",*indices))]:
        coefficients[name], pred = fitted(terms, ridge=5 if "indices" in name else 0)
        models[name] = pred
    def scores(group):
        y = np.array([z["B"] for z in group])
        return {name:float(np.sqrt(np.mean((y-pred(group))**2)))
                for name,pred in models.items()}
    return matched, {"source": BASE+"omni2_YEAR.dat",
                     "indices": ["Kp three-hourly (GFZ via OMNI)", "Dst hourly (Kyoto via OMNI)",
                                 "F10.7 daily (NRC via OMNI)"],
                     "daily_reducer": {"Kp": "UTC-day maximum", "Dst": "UTC-day minimum",
                                       "F10.7": "median of valid UTC-day hourly entries"},
                     "n_index_days": len(matched), "n_train": len(train), "n_test": len(test),
                     "correlations": r,
                     "storm_days": len(storm), "quiet_days": len(quiet),
                     "median_absolute_common_storm": (
                         float(np.median([z["abs_common"] for z in storm])) if storm else None),
                     "median_absolute_common_quiet": (
                         float(np.median([z["abs_common"] for z in quiet])) if quiet else None),
                     "train_coefficients": coefficients,
                     "train_rmse": scores(train), "heldout_rmse": scores(test),
                     "caveat": "Global indices are proxies, not local slant TEC or RM; observing schedule and solar-cycle trend are confounded."}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--directory", type=Path, default=HERE / ".omni_cache")
    ap.add_argument("--fetch", action="store_true",
                    help="Opt in to downloading NASA OMNI hourly files (2020–2023)")
    ap.add_argument("--write", action="store_true", help="Write diagnostic JSON and 124-day CSV")
    args = ap.parse_args()
    paths, shas = [], {}
    for year in range(2020, 2024):
        p, sha = fetch_year(year, args.directory, args.fetch)
        paths.append(p)
        shas[year] = sha
    matched, results = test_external_indices(indices_from_files(paths))
    results["sha256_omni2_year_files"] = shas
    if args.write:
        (HERE / "space_weather_checks.json").write_text(
            json.dumps(results, indent=2, sort_keys=True) + "\n")
        with (HERE / "space_weather_days.csv").open("w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=list(matched[0]))
            writer.writeheader()
            writer.writerows(matched)
    print(json.dumps({"observed_days": len(matched),
                      "correlations": results["correlations"],
                      "heldout_rmse": results["heldout_rmse"],
                      "storm_days": results["storm_days"],
                      "quiet_days": results["quiet_days"]}))


if __name__ == "__main__":
    main()
