"""Replay the lens search and independent echo/flux checks from the accompanying tables."""
import csv
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
TARGET = "J0826+7002"
LAGS = np.arange(-600, 605, 5)


def read(name):
    with (ROOT / name).open() as stream:
        return list(csv.DictReader(stream))


def curves(rows, timecol):
    """Group the tables into (image, band) arrays ordered by observation time."""
    grouped = {}
    for row in rows:
        values = [float(row[k]) for k in (timecol, "magnitude", "median_catalog_error")]
        grouped.setdefault((row["image"], row["band"]), []).append(values)
    return {key: np.array(sorted(vals)).T for key, vals in grouped.items()}


def smooth(obs, grid, width):
    times, mag, err = obs
    scaled = (grid[:, None] - times[None, :]) / width
    weights = np.where(abs(scaled) < 2, np.exp(-0.5 * scaled**2) / (err**2 + 0.035**2), 0)
    total = weights.sum(axis=1)
    return np.sum(weights * mag, axis=1) / np.where(total > 0, total, np.nan)


def corr(a, b):
    a, b = a - a.mean(), b - b.mean()
    return np.sum(a * b) / np.sqrt(np.sum(a**2) * np.sum(b**2))


def scan(leading, following, lags=LAGS, baseline=140):
    """Compare a finely smoothed curve to a delayed one only on supported grid points."""
    lo = min(leading[0].min(), following[0].min())
    hi = max(leading[0].max(), following[0].max())
    grid = 5 * np.arange(int(lo // 5) - 1, int(hi // 5) + 2)
    a, b = [smooth(x, grid, 10) for x in (leading, following)]
    if baseline:
        a -= smooth(leading, grid, baseline)
        b -= smooth(following, grid, baseline)
    result = []
    for lag in lags:
        later = np.roll(b, -int(lag // 5))
        valid = np.isfinite(a + later) & (grid + lag > grid[0]) & (grid + lag < grid[-1])
        result.append(corr(a[valid], later[valid]) if valid.sum() > 20 else np.nan)
    return np.array(result)


def print_peak(label, x, y, lags=LAGS, baseline=140):
    scores = scan(x, y, lags, baseline)
    if np.isfinite(scores).any():
        p = np.nanargmax(scores)
        at = lambda lag: scores[np.flatnonzero(lags == lag)[0]]
        print(f"{label:32s} {len(x[0]):4d} {len(y[0]):4d} {lags[p]:5d} "
              f"{scores[p]:7.3f} {at(210):7.3f} {at(0):7.3f}")
    return scores


def unique_pairs(t_b, t_a, lag, tolerance=12):
    distance = np.abs(t_a[:, None] - t_b[None, :] - lag)
    i_a, i_b = np.where(distance <= tolerance)
    ranked = sorted((distance[i, j], i, j) for i, j in zip(i_a, i_b))
    used_a, used_b, pairs = set(), set(), []
    for _, i, j in ranked:
        if i not in used_a and j not in used_b:
            used_a.add(i)
            used_b.add(j)
            pairs.append((j, i))
    return np.array(pairs, dtype=int).reshape(-1, 2).T


def red_noise(obs, tau, count, rng):
    """Independent unit-variance OU processes, sampled at the real observation times."""
    t, m, e = obs
    simulated = np.empty((count, len(t)))
    simulated[:, 0] = rng.normal(size=count)
    for i, gap in enumerate(np.diff(t), 1):
        rho = np.exp(-gap / tau)
        simulated[:, i] = rho * simulated[:, i - 1] + np.sqrt(1 - rho**2) * rng.normal(size=count)
    return np.std(m) * simulated + rng.normal(size=simulated.shape) * e


def red_noise_test(b, a):
    ib, ia = unique_pairs(b[0], a[0], 210)
    true = corr(b[1, ib], a[1, ia])
    for tau in (30, 100, 300, 1000):
        rng = np.random.default_rng(8130 + tau)
        x = red_noise(b, tau, 20000, rng)[:, ib]
        y = red_noise(a, tau, 20000, rng)[:, ia]
        x -= x.mean(1, keepdims=True)
        y -= y.mean(1, keepdims=True)
        correlations = (x * y).sum(1) / np.sqrt((x*x).sum(1) * (y*y).sum(1))
        exceeds = (correlations >= true).sum()
        print(f"fixed-lag null tau={tau:4} days: {exceeds}/20000 >= observed correlation")


def panstarrs():
    # Remove a free median for each band and image before pooling sparse noncontemporaneous bands.
    rows = read("ps1-detections.csv")
    data = {}
    for im in "BA":
        nights = {}
        for band in range(1, 6):
            selected = [r for r in rows if r["image"] == im and int(r["filter_id"]) == band]
            mag = np.array([8.9 - 2.5*np.log10(float(r["psf_jy"])) for r in selected])
            for r, m in zip(selected, mag - np.median(mag)):
                t = float(r["mjd"])
                err = 1.086 * float(r["error_jy"]) / float(r["psf_jy"])
                nights.setdefault(int(t), []).append([t, m, err])
        data[im] = np.array([np.median(nights[k], axis=0) for k in sorted(nights)]).T
    print_peak("PS1 pooled (exploratory)", data["B"], data["A"])
    for band in range(1, 6):
        times = {}
        for im in "BA":
            bins = {}
            for row in rows:
                if row["image"] == im and int(row["filter_id"]) == band:
                    t = float(row["mjd"])
                    bins.setdefault(int(t), []).append(t)
            times[im] = np.array([np.median(v) for _, v in sorted(bins.items())])
        ib, _ = unique_pairs(times["B"], times["A"], 210)
        print(f"PS1 filter {band}: {len(ib)} same-filter nightly pairs at +210 d (tolerance 12 d)")


def main():
    ztf, gaia = read("ztf-nights.csv"), curves(read("gaia-nights.csv"), "bjd_tcb_minus_2400000p5")
    print("run                              N(B) N(A) lag     max   r(210)    r(0)")
    for name in sorted({r["system"] for r in read("targets.csv")}):
        data = curves([r for r in ztf if r["system"] == name], "mjd")
        if not data:
            print(name, "no usable pairs")
        for band in ["zg", "zr"]:
            if ("A", band) in data and ("B", band) in data:
                print_peak(name + " " + band, data["B", band], data["A", band])

    data = curves([r for r in ztf if r["system"] == TARGET], "mjd")
    for band in ["zg", "zr", "zi"]:
        b, a = data["B", band], data["A", band]
        for t0, t1 in [(0, 59600), (59600, 1e6)]:
            early = lambda arr: arr[:, (arr[0] >= t0) & (arr[0] < t1)]
            print_peak(f"ZTF {band} {int(t0)}-{int(t1)}", early(b), early(a))
        for baseline in (0, 70, 140, 300):
            print_peak(f"ZTF {band} long scan base={baseline}", b, a, np.arange(-1200, 1205, 5), baseline)
    for band in ["g", "bp", "rp"]:
        b, a = gaia["B", band], gaia["A", band]
        for baseline in (0, 70, 140, 300):
            print_peak(f"Gaia {band} base={baseline}", b, a, baseline=baseline)
    b, a = gaia["B","g"], gaia["A","g"]
    for delay in (0, 210, -210, 155, 265, 365, -155, 575):
        ib, ia = unique_pairs(b[0], a[0], delay)
        print(f"Gaia raw G at {delay:+d} d: {len(ib)} one-to-one pairs, r={corr(b[1,ib],a[1,ia]):.3f}")
    for block in range(5):
        omit = lambda x: x[:, (x[0] - a[0].min()) // 220 != block]
        scores = scan(omit(b), omit(a))
        print(f"Gaia G omit 220-day block {block + 1}: best lag {LAGS[np.nanargmax(scores)]} d")
    red_noise_test(b, a)
    panstarrs()
    for band in ["zg","zr","zi"]:
        rows = [r for r in read("field-controls.csv") if r["band"] == band]
        changes = {r["object"]: float(r["late_minus_early_mag"]) for r in rows if r["object"] != "field"}
        field = [abs(float(r["late_minus_early_mag"])) for r in rows if r["object"]=="field"]
        print(f"{band} A {changes['A']:+.3f}, B {changes['B']:+.3f}, A-B {changes['A']-changes['B']:+.3f}, "
              f"{len(field)} field comparisons, max absolute change {max(field):.3f}")
    alerts = read("difference-detections.csv")
    for lo, hi in [(58900,59500),(60350,60970)]:
        selected = [r for r in alerts if lo <= float(r["mjd"]) < hi]
        npos = sum(float(r["isdiffpos"])>0 for r in selected)
        print(f"MJD [{lo},{hi}): {npos}/{len(selected)} positive difference detections, median separation "
              f"{np.median([float(r['off']) for r in selected]):.3f} arcsec")


if __name__ == "__main__":
    main()
