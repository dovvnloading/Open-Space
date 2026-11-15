#!/usr/bin/env python3
"""New source/epoch holdouts and ionosphere-origin tests for Keith et al. MeerKAT RMs.

No network required. Reads PR #47 inputs without changing them.
Run: python results/meerkat-common-mode-origin/diagnostics.py --write
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / "meerkat-shared-faraday-rotation"
MEASUREMENT_SHA256 = "0536aac6f163dfba2191f8728a073a2fb9b38c03ce33425376091172cbd490ab"
LAT, LON = -30.71322, 21.44306


def prior_replay():
    """Import the frozen PR #47 selection, not an independently altered pipeline."""
    spec = importlib.util.spec_from_file_location("meerkat_pr47_replay", PRIOR / "replay.py")
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def dates_mjd_year(mjd):
    return (dt.date(1858, 11, 17) + dt.timedelta(days=int(mjd))).year


def read_original_csv(name):
    with (PRIOR / name).open(newline="") as f:
        return list(csv.DictReader(f))


def matched_peers(prior, stars):
    """Exactly the original angular/time/precision selection, retaining all cases."""
    sessions = prior.peer_table(stars)
    min_dot = np.cos(np.deg2rad(30.0))
    cases = {}
    for name, (t, y, e, sky, ion, raw) in stars.items():
        entries = []
        for i, stamp in enumerate(t):
            if e[i] >= 2:
                continue
            peers = defaultdict(list)
            for neighbour, tp, dep, direction in sessions[int(stamp)]:
                if (neighbour != name and abs(tp - stamp) < 2.16 / 24
                        and sky @ direction > min_dot):
                    peers[neighbour].append(dep)
            if len(peers) >= 5:
                prediction = float(np.median([np.median(v) for v in peers.values()]))
                entries.append((stamp, y[i], e[i], prediction, ion[i], raw[i], len(peers)))
        if (len(entries) >= 18 and
                entries[-1][0] - entries[0][0] >= 800):
            cases[name] = np.asarray(entries)
    return cases


def ionosphere_trends(prior, data):
    time, corrected, err, peers, ion, raw, _ = data.T
    slopes = {}
    for key, y in [("observed_rm", raw),
                   ("published_iono", ion),
                   ("post_iono_rm", corrected),
                   ("peer_prediction", peers),
                   ("post_peer_rm", corrected - peers)]:
        slope, se, _ = prior.line(time, y, err)
        slopes[key] = {"slope": float(slope), "conditional_se": float(se)}
    return {"n": len(data), "slopes": slopes,
            "peer_iono_r": pearson(peers, ion),
            "peer_iono_beta": float(np.cov(peers, ion, bias=True)[0, 1] / np.var(ion))}


def sky_sector(cases, stars):
    per_source = {}
    for name, data in cases.items():
        time, _, _, p, ion, _, _ = data.T
        sky = stars[name][3]
        ra_h = (np.arctan2(sky[1], sky[0]) % (2 * np.pi)) * 12 / np.pi
        dec = np.rad2deg(np.arcsin(sky[2]))
        centered_x, centered_y = ion - np.median(ion), p - np.median(p)
        beta = float(centered_x @ centered_y / (centered_x @ centered_x))
        early, late = time < 59580, time >= 59580
        per_source[name] = {
            "n": len(time), "ra_h": float(ra_h), "dec_deg": float(dec),
            "r_peer_iono": pearson(p, ion), "beta_peer_iono": beta,
            "r_early": pearson(p[early], ion[early]) if early.sum() >= 7 else None,
            "r_late": pearson(p[late], ion[late]) if late.sum() >= 7 else None}
    bins = {}
    for lower in [3, 6, 9, 12, 15, 18, 21]:
        group = [z for z in per_source.values() if lower <= z["ra_h"] < lower + 3]
        if not group:
            continue
        bins[f"{lower:02d}-{lower+3:02d}h"] = {
            "n": len(group),
            "median_r": float(np.median([z["r_peer_iono"] for z in group])),
            "median_beta": float(np.median([z["beta_peer_iono"] for z in group])),
            "r_below_minus_half": sum(z["r_peer_iono"] < -0.5 for z in group),
            "r_above_plus_half": sum(z["r_peer_iono"] > 0.5 for z in group)}
    paired = [z for z in per_source.values()
              if z["r_early"] is not None and z["r_late"] is not None]
    summary = {
        "matched_sources": len(per_source),
        "n_r_below_minus_half": sum(z["r_peer_iono"] < -0.5 for z in per_source.values()),
        "n_r_above_plus_half": sum(z["r_peer_iono"] > 0.5 for z in per_source.values()),
        "early_late_enough": len(paired),
        "early_late_sign_concordant": sum(
            np.sign(z["r_early"]) == np.sign(z["r_late"]) for z in paired),
        "ra_bins": bins}
    return per_source, summary


def solar_altitude(mjd, ra_deg, dec_deg):
    """Approximate UTC mean solar hour, topocentric geometric elevation in degrees."""
    lst = (18.697374558 + 24.06570982441908 * (mjd - 51544.5)
           + LON / 15) % 24
    ha_hours = ((lst - ra_deg / 15 + 12) % 24) - 12
    h = np.deg2rad(ha_hours * 15)
    lat, dec = np.deg2rad([LAT, dec_deg])
    elevation = np.rad2deg(np.arcsin(
        np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(h)))
    solar_hour = (mjd % 1 * 24 + LON / 15) % 24
    return elevation, solar_hour, ha_hours


def precompute_observations(stars):
    all_rows = []
    for name, (time, y, e, sky, ion, raw) in stars.items():
        ra = (np.arctan2(sky[1], sky[0]) % (2 * np.pi)) * 180 / np.pi
        dec = np.rad2deg(np.arcsin(sky[2]))
        for t, residual, sigma, z in zip(time, y, e, ion):
            if sigma >= 1.8:
                continue
            elevation, solar_hour, ha = solar_altitude(t, ra, dec)
            all_rows.append((name, int(t), residual, z,
                             elevation, solar_hour, ha))
    # Each observation is one pulsar-epoch, not a distinct independent test.
    by_day = defaultdict(set)
    for name, day, *_ in all_rows:
        by_day[day].add(name)
    days = {d for d, names in by_day.items() if len(names) >= 12}
    return [v for v in all_rows if v[1] in days]


def fixed_effect_residual(y, source_ids, day_ids):
    """Alternating projection onto orthogonal complement of source/day dummies."""
    r = np.asarray(y, float).copy()
    for _ in range(50):
        step = 0.0
        for ids in (source_ids, day_ids):
            sums = np.bincount(ids, weights=r)
            count = np.bincount(ids)
            shift = (sums / count)[ids]
            r -= shift
            step = max(step, float(np.max(np.abs(shift))))
        if step < 1e-8:
            break
    return r


def fixed_effect_model(rows, features=("iono",)):
    """Source+date FE with date-clustered sandwich SE: *association*, not causal slope."""
    by_day = defaultdict(set)
    for v in rows:
        by_day[v[1]].add(v[0])
    kept = [v for v in rows if len(by_day[v[1]]) >= 12]
    if len(kept) < 200:
        return None
    source = np.unique([v[0] for v in kept], return_inverse=True)[1]
    day = np.unique([v[1] for v in kept], return_inverse=True)[1]
    residual = np.array([v[2] for v in kept])
    iono = np.array([v[3] for v in kept])
    alt = np.array([v[4] for v in kept])
    sun = np.array([v[5] for v in kept])
    ha = np.array([v[6] for v in kept])
    value = {"iono": iono, "sec": 1 / np.sin(np.deg2rad(alt)),
             "sinSun": np.sin(sun * np.pi / 12),
             "cosSun": np.cos(sun * np.pi / 12),
             "sinHA": np.sin(ha * np.pi / 12),
             "cosHA": np.cos(ha * np.pi / 12)}
    y = fixed_effect_residual(residual, source, day)
    X = np.column_stack([fixed_effect_residual(value[k], source, day)
                         for k in features])
    xx = X.T @ X
    inv = np.linalg.inv(xx)
    beta = inv @ (X.T @ y)
    e = y - X @ beta
    scores = np.zeros((day.max() + 1, len(features)))
    np.add.at(scores, day, X * e[:, None])
    cov = inv @ (scores.T @ scores) @ inv * (len(scores) / (len(scores) - 1))
    return {"n_visits": len(kept), "n_pulsars": int(source.max()+1),
            "n_dates": len(scores),
            "terms": {name: {"beta": float(beta[i]),
                             "day_cluster_se": float(np.sqrt(max(0, cov[i,i])))}
                      for i, name in enumerate(features)},
            "partial_R2": float(1 - e @ e / (y @ y))}


def heldout_calendar_test():
    """A and B are distinct SHA-256 pulsar sets from the frozen prior pipeline."""
    dates = read_original_csv("dates.csv")
    r = []
    for z in dates:
        year = dates_mjd_year(int(float(z["mjd_utc"])))
        r.append(dict(year=year, a=float(z["resid_A"]),
                      b=float(z["resid_B"]), ion=float(z["iono_A"])))
    train, test = [z for z in r if z["year"] <= 2021], [z for z in r if z["year"] >= 2022]
    assert len(train) == 56 and len(test) == 68
    def fit(terms):
        X = np.array([[1, *[z[k] for k in terms]] for z in train])
        y = np.array([z["b"] for z in train])
        return np.linalg.lstsq(X, y, rcond=None)[0]
    coefficients = {"peers": fit(["a"]), "ionosphere": fit(["ion"]),
                    "peers_and_ionosphere": fit(["a", "ion"])}
    def errors(samples):
        y = np.array([z["b"] for z in samples])
        a = np.array([z["a"] for z in samples])
        ion = np.array([z["ion"] for z in samples])
        predictions = {"no_correction": np.zeros(len(y)),
                       "unscaled_peers": a,
                       "trained_peers": coefficients["peers"][0] + coefficients["peers"][1] * a,
                       "trained_ionosphere": coefficients["ionosphere"][0] + coefficients["ionosphere"][1] * ion,
                       "trained_peers_plus_iono": (coefficients["peers_and_ionosphere"][0]
                           + coefficients["peers_and_ionosphere"][1] * a
                           + coefficients["peers_and_ionosphere"][2] * ion)}
        return {k: {"rmse": float(np.sqrt(np.mean((y - p)**2))),
                    "mae": float(np.mean(np.abs(y - p)))}
                for k, p in predictions.items()}
    scores_train, scores_test = errors(train), errors(test)
    # Fixed pseudorandom sequence, resampling *dates*, not thousands of correlated visits.
    state, gains = 81321, []
    for _ in range(10_000):
        ss0 = ss1 = 0.
        for _ in test:
            state = (1664525 * state + 1013904223) & 0xffffffff
            z = test[int((state / 2**32) * len(test))]
            ss0 += z["b"]**2
            ss1 += (z["b"] - z["a"])**2
        gains.append(1 - np.sqrt(ss1 / ss0))
    gains.sort()
    return {"n_training_days": len(train), "n_heldout_days": len(test),
            "train_period": "2020-2021", "test_period": "2022-2023",
            "train_beta": {k: [float(q) for q in v] for k,v in coefficients.items()},
            "train": scores_train, "heldout": scores_test,
            "heldout_peer_correlation": pearson([z["a"] for z in test], [z["b"] for z in test]),
            "heldout_uncalibrated_rmse_gain": 1 - (scores_test["unscaled_peers"]["rmse"]
                                                  / scores_test["no_correction"]["rmse"]),
            "heldout_day_bootstrap_95pct": [gains[250], gains[9749]]}


def published_risk_rows(info):
    """Audit the full-series >3σ selection without claiming all 76 are covered."""
    slopes = read_original_csv("slopes.csv")
    result = []
    for row in slopes:
        d = {k: float(v) for k,v in row.items() if k != "psr"}
        source = row["psr"]
        if source not in info:
            raise AssertionError(f"Matched pulsar missing from new peer diagnostic: {source}")
        full = abs(d["full_slope"]) > 3 * d["full_se"]
        before = abs(d["matched_slope"]) > 3 * d["matched_se"]
        after = abs(d["differential_slope"]) > 3 * d["differential_se"]
        if full and not after:
            label = "full_series_3sigma_lost"
        elif full and after:
            label = "full_series_3sigma_persists"
        elif before and not after:
            label = "matched_series_3sigma_lost"
        elif not before and after:
            label = "matched_series_3sigma_gained"
        else:
            label = "not_nominally_3sigma_in_scope"
        result.append(dict(psr=source, n_matched=int(d["n_matched"]),
                           full_slope=d["full_slope"], full_se=d["full_se"],
                           matched_slope=d["matched_slope"], matched_se=d["matched_se"],
                           differential_slope=d["differential_slope"],
                           differential_se=d["differential_se"],
                           removed_slope=d["matched_slope"] - d["differential_slope"],
                           full_3sigma=int(full), matched_3sigma=int(before),
                           differential_3sigma=int(after), risk_class=label,
                           ra_h=info[source]["ra_h"], dec_deg=info[source]["dec_deg"],
                           peer_iono_r=info[source]["r_peer_iono"],
                           peer_iono_beta=info[source]["beta_peer_iono"]))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true",
                        help="Write new-directory checks.json and published_risks.csv")
    args = parser.parse_args()
    actual = hashlib.sha256((PRIOR / "measurements.csv").read_bytes()).hexdigest()
    if actual != MEASUREMENT_SHA256:
        raise ValueError(f"Unexpected input checksum: {actual}")
    prior = prior_replay()
    stars = prior.read_rows()
    cases = matched_peers(prior, stars)
    assert len(cases) == 304 and len(cases["J1057-5226"]) == 53
    source, sky = sky_sector(cases, stars)
    obs = precompute_observations(stars)
    fe = {"all_iono": fixed_effect_model(obs, ("iono",)),
          "all_iono_and_geometry": fixed_effect_model(obs, ("iono", "sec", "sinSun", "cosSun")),
          "daylight_iono": fixed_effect_model([z for z in obs if 8 <= z[5] < 16], ("iono",)),
          "other_hours_iono": fixed_effect_model([z for z in obs if not (8 <= z[5] < 16)], ("iono",))}
    holdout = heldout_calendar_test()
    risks = published_risk_rows(source)
    report = {"provenance": {"input_sha256": actual, "source_rows": 20_365,
                             "matching": "frozen PR47 selection; target excluded from peers"},
              "J1057": ionosphere_trends(prior, cases["J1057-5226"]),
              "sky_sector": sky, "fixed_effects": fe,
              "year_holdout": holdout,
              "published_risk_scope": {
                  "matched_sources": len(risks),
                  "full_series_3sigma": sum(z["full_3sigma"] for z in risks),
                  "full_series_3sigma_to_subthreshold": sum(
                      z["full_3sigma"] and not z["differential_3sigma"] for z in risks),
                  "matched_series_3sigma": sum(z["matched_3sigma"] for z in risks),
                  "matched_series_3sigma_to_subthreshold": sum(
                      z["matched_3sigma"] and not z["differential_3sigma"] for z in risks),
                  "new_matched_3sigma": sum(
                      not z["matched_3sigma"] and z["differential_3sigma"] for z in risks)}}
    if args.write:
        (HERE / "checks.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=lambda z: z.item() if isinstance(z, np.generic) else str(z)) + "\n")
        with (HERE / "published_risks.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(risks[0]))
            writer.writeheader()
            writer.writerows(risks)
    # Compact console summary; full results saved on --write.
    print(json.dumps({
        "J1057_slopes": {k: round(v["slope"], 6)
                         for k, v in report["J1057"]["slopes"].items()},
        "sky_sector": {k: sky[k] for k in
                       ["matched_sources", "n_r_below_minus_half", "n_r_above_plus_half"]},
        "ionosphere_beta_FE": round(fe["all_iono"]["terms"]["iono"]["beta"], 4),
        "heldout_dates": holdout["n_heldout_days"],
        "heldout_rmse": {k: round(v["rmse"], 5)
                         for k,v in holdout["heldout"].items()},
        "published_risk_scope": report["published_risk_scope"]}, default=lambda z: z.item() if isinstance(z, np.generic) else str(z)))


if __name__ == "__main__":
    main()
