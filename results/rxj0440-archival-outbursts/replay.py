"""Reproduce the event search and the GBM-selected tests from the frozen tables."""
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OFF = np.r_[np.ones(20), np.zeros(11), np.ones(20)]


def load(name):
    return np.loadtxt(DATA / (name + "-daily.csv.gz"), delimiter=",", skiprows=1)


def daily_scores(rows, days):
    """Compare three consecutive days to days 6--25 on both sides."""
    weights = np.zeros(len(days))
    signal = np.zeros(len(days))
    counts = np.zeros(len(days))
    indices = np.floor(rows[:, 0]).astype(int) - days[0]
    have = (indices >= 0) & (indices < len(days))
    idx = indices[have]
    if len(idx) != len(set(idx)):
        raise ValueError("Duplicate day labels")
    weights[idx] = 1 / rows[have, 2] ** 2
    signal[idx] = rows[have, 1] * weights[idx]
    counts[idx] = 1
    block = np.ones(3)
    w = np.convolve(weights, block, "same")
    b = np.convolve(weights, OFF, "same")
    available = ((w > 0) & (b > 0) & (np.convolve(counts, OFF, "same") >= 8))
    z = np.full(len(days), np.nan)
    q = np.convolve(signal, block, "same")[available] / w[available]
    base = np.convolve(signal, OFF, "same")[available] / b[available]
    z[available] = (q - base) / np.sqrt(1 / w[available] + 1 / b[available])
    return z


def episodes(days, zb, zm):
    passes = np.flatnonzero((zb > 4) & (zm > 4))
    if not len(passes):
        return []
    groups = np.split(passes, np.flatnonzero(np.diff(passes) > 12) + 1)
    return [int(g[np.argmax(np.minimum(zb[g], zm[g]))]) for g in groups]


def weighted(rows):
    if not len(rows):
        return None
    w = rows[:, 2] ** -2
    return float(w @ rows[:, 1] / w.sum()), float(1 / np.sqrt(w.sum()))


def contrast(rows, start, stop):
    t = rows[:, 0]
    center = (start + stop) / 2
    during = rows[(t >= start) & (t < stop)]
    baseline = rows[(abs(t - center) < 30) & (abs(t - center) >= 5)]
    if not len(during) or not len(baseline):
        return None
    y, e = weighted(during)
    y0, e0 = weighted(baseline)
    peak = during[np.argmax(during[:, 1])]
    return dict(n=len(during), baseline_n=len(baseline), mean=y, error=e,
                baseline=y0, baseline_error=e0, z=(y - y0) / np.hypot(e, e0),
                brightest_day_label=float(peak[0]), brightest_day_rate=float(peak[1]))


def run():
    bat = load("bat")
    maxi = load("maxi")
    bat = bat[(bat[:, 3] == 0) & (bat[:, 2] > 0)]
    maxi = maxi[(maxi[:, 2] > 0) & (maxi[:, 2] < 0.2)]
    gbm = np.genfromtxt(DATA / "gbm-detections.csv", names=True, delimiter=",")
    gbm = gbm[gbm["detected"] == 1]
    joint_start = int(max(bat[0, 0], maxi[0, 0]))
    joint_stop = int(min(bat[-1, 0], maxi[-1, 0]))
    days = np.arange(joint_start, joint_stop + 1)
    zb, zm = daily_scores(bat, days), daily_scores(maxi, days)
    peaks = episodes(days, zb, zm)
    events = [dict(mjd=int(days[i]), bat_z=float(zb[i]), maxi_z=float(zm[i])) for i in peaks]
    # Shift the complete MAXI statistic with its own window, keeping BAT fixed.
    # Known major flares remain in both series. This is a coincidence diagnostic.
    shifts = np.arange(30, len(days) - 29)
    null = np.array([len(episodes(days, zb, np.roll(zm, int(s)))) for s in shifts])
    first = gbm[0]
    recent = gbm[-1]
    sep = gbm[(gbm["barytime_mjd"] >= 61300) & (gbm["barytime_mjd"] < 61330)][0]
    targeted = {}
    for name, r in [("february_2010", first), ("november_2026", recent)]:
        f, e = r["frequency_hz"], r["frequency_error_hz"]
        lo, hi = r["range_start_mjd"], r["range_stop_mjd"]
        targeted[name] = dict(range_mjd=[float(lo),float(hi)], epoch_mjd=float(r["barytime_mjd"]),
                              period_s=1/f, period_error_s=e/f**2,
                              amplitude_over_error=r["amplitude"]/r["amplitude_error"],
                              bat=contrast(bat,lo,hi), maxi=contrast(maxi,lo,hi))
    return dict(joint_scan_day_labels=[joint_start,joint_stop],
                available_simultaneous_centers=int(np.isfinite(zb+zm).sum()),
                episodes=events, shift_count=len(null), median_shifted_count=float(np.median(null)),
                maximum_shifted_count=int(null.max()), shifts_matching_observed=int((null>=len(peaks)).sum()),
                fixed_gbm_windows=targeted,
                february_to_april_first_gbm_epochs_d=float(gbm[1]["barytime_mjd"]-first["barytime_mjd"]),
                september_to_november_first_gbm_epochs_d=float(recent["barytime_mjd"]-sep["barytime_mjd"]))


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
