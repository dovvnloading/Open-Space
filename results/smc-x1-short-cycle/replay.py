"""Measure the cycle lengths, timing ranks and independent checks."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parent
# Windows were chosen during exploration. These are not blinded detections.
WINDOWS = {"before": (60490, 60695), "fast": (60695, 60860),
           "after": (60865, 61080), "2026": (61136, 61349)}


def data():
    path = ROOT / "data/measurements.csv"
    manifest = json.loads((ROOT / "data/sources.json").read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest() != manifest["measurements_sha256"]:
        raise ValueError("Measurement hash differs from sources.json")
    result = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            result.setdefault(row["series"], []).append([float(row[k]) for k in
                              ("mjd", "rate", "error", "n_measurements")])
    return {key: np.asarray(values) for key, values in result.items()}


def design(t, period, harmonics=1):
    angle = 2 * np.pi * (t - 60700) / period
    return np.column_stack([np.ones(len(t))] + [f(h * angle)
                           for h in range(1, harmonics + 1) for f in (np.cos, np.sin)])


def fit(a, period, harmonics=1):
    t, y, err = a[:, :3].T
    x = design(t, period, harmonics)
    b = np.linalg.lstsq(x / err[:, None], y / err, rcond=None)[0]
    return np.sum(((y - x @ b) / err)**2), b


def optimum(a, harmonics=1):
    # The common, wide grid also allows the 30-day control as a blind peak.
    grid = np.arange(24, 76.01, .25)
    i = np.argmin([fit(a, p, harmonics)[0] for p in grid])
    if i in (0, len(grid) - 1):
        return float(grid[i])
    return minimize_scalar(lambda p: fit(a, p, harmonics)[0],
                           bounds=(grid[i-1], grid[i+1]), method="bounded").x


def select(a, lo, hi):
    return a[(a[:, 0] >= lo) & (a[:, 0] < hi)]


def cycle_minima(a, width, prominence):
    day = a[:, 0].astype(int)
    origin = day.min()
    weights = 1 / a[:, 2]**2
    norm = gaussian_filter1d(np.bincount(day-origin, weights), float(width))
    flux = gaussian_filter1d(np.bincount(day-origin, weights * a[:, 1]), float(width))
    model = np.divide(flux, norm, out=np.zeros_like(norm), where=norm > 0)
    indices, _ = find_peaks(-model, distance=29, prominence=prominence)
    return indices + origin + .5


def rank_runs(a, width):
    times = cycle_minima(a, width, .0023)
    eligible = []
    for i in range(len(times) - 4):
        dates = times[i:i+5]
        if np.diff(dates).max() > 75:
            continue
        obs = select(a, dates[0] - 3, dates[-1] + 3.01)[:, 0]
        if not len(obs) or max(np.min(abs(obs - t)) for t in dates) > 5:
            continue
        if np.diff(np.r_[dates[0], obs, dates[-1]]).max() > 12:
            continue
        eligible.append((float((dates[-1] - dates[0]) / 4), dates.tolist()))
    eligible.sort()
    return eligible


def segment(a, lo, hi):
    b = select(a, lo, hi)
    n = len(b)
    if n < 20:
        return {"days": n}
    largest_gap = np.diff(np.r_[lo, b[:, 0], hi]).max()
    result = {"days": n, "max_gap_d": round(float(largest_gap), 1)}
    if largest_gap <= .25 * (hi-lo):
        result["period_d"] = round(float(optimum(b)), 3)
    return result


def evaluate_template(reference, target):
    p = optimum(reference)
    _, b = fit(reference, p)
    t, y, e = target[:, :3].T
    h = design(t, p)[:, 1:] @ b[1:]
    def score(x):
        coef = np.linalg.lstsq(x / e[:, None], y / e, rcond=None)[0]
        r = y - x @ coef
        return float(coef[-1]), float(np.sqrt(np.average(r*r, weights=e**-2)))
    return {"period_d": round(float(p), 3), "gain_and_weighted_rms":
            [round(v, 6) for v in score(np.column_stack([np.ones(len(t)), h]))],
            "constant_weighted_rms": round(score(np.ones((len(t), 1)))[1], 6)}


def projected_minimum(a, lo, harmonics):
    b = select(a, lo, 61349)
    p = optimum(b, harmonics)
    _, coefficients = fit(b, p, harmonics)
    days = np.arange(61358., 61395., .05)
    curve = design(days, p, harmonics) @ coefficients
    return round(float(days[np.argmin(curve)]), 2)


def measure():
    series = data()
    table = {key: {epoch: segment(a, *window) for epoch, window in WINDOWS.items()}
             for key, a in series.items() if not key.endswith("_core") and key != "smcfermi"}
    rankings = {}
    for width in range(2, 7):
        ordered = rank_runs(series["smcbat"], width)
        matches = [(i + 1, d, dates) for i, (d, dates) in enumerate(ordered)
                   if 60670 < dates[0] < 60725 and dates[-1] < 60865]
        rankings[width] = {"eligible_runs": len(ordered), "three_fastest": ordered[:3],
                           "event_rank_period_minima": matches}
    swift = select(series["smcbat"], *WINDOWS["fast"])
    maxi = select(series["smcmaxi"], *WINDOWS["fast"])
    predictions = {"Swift_to_MAXI": evaluate_template(swift, maxi)}
    gbm = select(series["smcfermi"], *WINDOWS["fast"])
    predictions["Swift_to_GBM_preliminary"] = {
        **evaluate_template(swift, gbm), "four_day_bins": len(gbm),
        "free_period_d": round(float(optimum(gbm)), 3)}
    for region in ("outside", "core"):
        target = select(series["smc_maxi_orb_" + region], *WINDOWS["fast"])
        predictions["Swift_to_MAXI_" + region] = evaluate_template(swift, target)
    mask_checks = {}
    for k in ("smcbat", "smcmaxi", "smc_bat_orb_outside", "smc_maxi_orb_outside"):
        a = series[k]
        b = select(a, *WINDOWS["fast"])
        periods = [optimum(b[((b[:, 0]-60695)//38.5).astype(int) != j]) for j in range(4)]
        mask_checks[k] = {"two_harmonics_d": round(float(optimum(b, 2)), 3),
                          "omit_each_of_four_38p5d_blocks_d": [round(float(v), 3) for v in periods],
                          "shift_both_window_edges_d": {
                              shift: round(float(optimum(select(a, 60695+shift, 60860+shift))), 3)
                              for shift in (-10, -5, 5, 10)}}
    minima = {}
    for m, key, prominence in (("Swift", "smcbat", .0023), ("MAXI", "smcmaxi", .04)):
        dates = cycle_minima(series[key], 3, prominence)
        local = dates[(dates >= 60690) & (dates < 60860)]
        obs = select(series[key], 60690, 60860)[:,0]
        minima[m] = {"mjd": local.tolist(), "intervals_d": np.diff(local).tolist(),
                     "distance_nearest_observation_d": [float(np.min(abs(obs - t))) for t in local]}
    forecast = {f"{lo}_h{h}": projected_minimum(series["smcmaxi"], lo, h)
                for lo in (61136, 61180, 61210) for h in (1, 2)}
    archive_runs = rank_runs(series["smcbat"], 3)
    return dict(windows_mjd=WINDOWS, timing=table, minima=minima,
                projected_minimum_mjd=forecast, swift_archive_rank=rankings,
                eligible_three_day_smoothed_runs=archive_runs,
                robustness=mask_checks, cross_mission_waveform=predictions)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = json.dumps(measure(), indent=2) + "\n"
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end="")
