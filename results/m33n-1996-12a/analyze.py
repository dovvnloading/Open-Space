"""Reproduce the Williams & Shafter position matches and modern photometry."""
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
ARCSEC = 180.0 * 3600 / np.pi


def read(name):
    with (ROOT / name).open() as stream:
        return list(csv.DictReader(stream))


def units(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.column_stack([np.cos(ra) * np.cos(dec), np.sin(ra) * np.cos(dec), np.sin(dec)])


def angles(chord):
    return 2 * ARCSEC * np.arcsin(np.clip(chord / 2, 0, 1))


def spatial():
    positions = read("original_positions.csv")
    variables = read("mira_fields.csv")
    ra = np.array([float(p["ra_deg"]) for p in positions])
    dec = np.array([float(p["dec_deg"]) for p in positions])
    rng = np.random.default_rng(87544)
    repetitions = 100000
    azimuth = rng.uniform(0, 2 * np.pi, repetitions)
    radius = np.sqrt(rng.uniform(12 ** 2, 60 ** 2, repetitions))
    # A common local east/north shift of all six pointings preserves their geometry.
    random_ra = ra[:, None] + radius * np.cos(azimuth) / (3600 * np.cos(np.radians(dec[:, None])))
    random_dec = dec[:, None] + radius * np.sin(azimuth) / 3600
    vectors = units(random_ra.ravel(), random_dec.ravel())
    result = {}
    for key, samples in [("all", variables), ("gold", [s for s in variables if s["classifier_votes"] == "6"])]:
        tree = cKDTree(units([float(s["ra_deg"]) for s in samples], [float(s["dec_deg"]) for s in samples]))
        ds, nearest = tree.query(units(ra, dec))
        ds = angles(ds)
        # A 75" data extraction margin suffices for all shifted searches for
        # matches within the 1.54" comparison threshold (maximum shift 60").
        scrambled = np.sort(angles(tree.query(vectors)[0]).reshape(len(ra), -1), axis=0)
        result[key] = {
            "objects_within_75arcsec": len(samples),
            "nearest": {p["id"]: [samples[i]["id"], round(float(d), 4)] for p, i, d in zip(positions, nearest, ds)},
            "offset_draws": repetitions,
            "any_one_at_least_as_close_as_1996": int(np.count_nonzero(scrambled[0] <= ds[3])),
            "any_two_as_close_as_1996_and_1997": int(np.count_nonzero((scrambled[0] <= ds[3]) & (scrambled[1] <= ds[4]))),
        }
    return result


def light_curves():
    data = read("ztf_dr24.csv")
    alerts = read("ztf_alerts.csv")
    output = {}
    for case in ["96", "97"]:
        points = [r for r in data if r["case"] == case and r["filtercode"] == "zi" and r["catflags"] == "0"]
        mags = [float(r["mag"]) for r in points]
        times = [float(r["mjd"]) for r in points]
        active = sorted(float(a["mjd"]) for a in alerts if a["case"] == case and a["isdiffpos"] == "1" and a["fid"] == "2")
        seasons = []
        for mjd in active:
            if not seasons or mjd - seasons[-1][-1] > 90:
                seasons.append([])
            seasons[-1].append(mjd)
        output[case] = {
            "ztf_i_catflags_zero": len(mags),
            "i_range_mag": [round(min(mags), 4), round(max(mags), 4)],
            "i_span_days": round(max(times) - min(times), 2),
            "positive_r_alert_groups": [[len(s), round(s[-1] - s[0], 2)] for s in seasons],
        }
    return output


def phase_sensitivity():
    # Compare surveys in i, allowing different zero points. Nightly aggregation
    # prevents four closely spaced Pan-STARRS exposures from dominating the fit.
    nightly = defaultdict(list)
    for r in read("ps1_detections.csv"):
        if r["filterID"] == "3":
            flux, error = float(r["psfFlux"]), float(r["psfFluxErr"])
            mjd = float(r["obsTime"])
            nightly[0, int(mjd)].append([mjd, -2.5 * np.log10(flux / 3631), 2.5 * error / (np.log(10) * flux)])
    for r in read("ztf_dr24.csv"):
        if r["case"] == "96" and r["filtercode"] == "zi" and r["catflags"] == "0":
            mjd = float(r["mjd"])
            nightly[1, int(mjd)].append([mjd, float(r["mag"]), float(r["magerr"])])
    t, m, err = np.array([np.median(r, axis=0) for r in nightly.values()]).T
    group = np.array([k[0] for k in nightly])
    weights = 1 / np.hypot(err, 0.12)  # model floor in mag, descriptive fits only

    def solve(period):
        angle = 2 * np.pi * (t - 59000) / period
        matrix = np.column_stack([group == 0, group == 1, np.cos(angle), np.sin(angle)])
        beta = np.linalg.lstsq(matrix * weights[:, None], m * weights, rcond=None)[0]
        phase = np.arctan2(-beta[-1], -beta[-2]) / (2 * np.pi)
        historical_phase = ((50426 - 59000) / period - phase + 0.5) % 1 - 0.5
        return float(np.sum(((m - matrix @ beta) * weights) ** 2)), float(historical_phase)

    def best(periods):
        return min((solve(p)[0], p, solve(p)[1]) for p in periods)

    fixed = 543.85
    short = best(np.arange(200, 351, 0.1))
    long = best(np.arange(450, 651, 0.1))
    return {
        "ps1_nights": int(np.sum(group == 0)),
        "ztf_nights": int(np.sum(group == 1)),
        "period_and_phase_relative_to_i_max": [
            [fixed, round(solve(fixed)[1], 3)],
            [round(float(long[1]), 1), round(long[2], 3)],
            [round(float(short[1]), 1), round(short[2], 3)],
        ],
    }


if __name__ == "__main__":
    results = {"positions": spatial(), "modern_photometry": light_curves(), "constant_period_sensitivity": phase_sensitivity()}
    print(json.dumps(results, indent=2))
