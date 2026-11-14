"""Measure the eclipse clock of each binary, with both apertures and time baselines."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
CLOCKS = {
    "A primary": (1.361339, 1440.2465, .016),
    "A secondary": (1.361339, 1440.2465 + 1.361339 / 2, .017),
    "B primary": (1.485117, 1443.958, .016),
}
SECTORS = (5, 32, 98, 109)


def one_clock(time, flux, name, half_window, cadence):
    period, epoch, width = CLOCKS[name]
    other_eclipse = np.zeros(len(time), dtype=bool)
    for label, (p, e, _) in CLOCKS.items():
        if label != name:
            phase = ((time - e) / p + .5) % 1 - .5
            other_eclipse |= abs(phase) <= .064
    cycles = np.arange(np.floor((time.min()-epoch)/period),
                       np.ceil((time.max()-epoch)/period)+1, dtype=int)
    for cycle in cycles:
        mid = epoch + period * cycle
        good = (abs(time - mid) < half_window) & ~other_eclipse
        x, y = time[good] - mid, flux[good]
        if sum(abs(x) < .03) < 2 or sum(x < -.05) < 2 or sum(x > .05) < 2:
            continue

        def residual(params):
            intercept, slope, depth, shift = params
            # Midpoint integration of the local Gaussian profile. A separate
            # baseline slope and eclipse depth absorb slow flux/shape changes.
            inside = x[:, None] + cadence * np.array([-.4, -.2, 0, .2, .4]) - shift
            profile = np.exp(-.5 * (inside / (width * period))**2).mean(axis=1)
            return intercept + slope*x - depth*profile - y

        fit = least_squares(residual, [1, 0, .08, 0],
                            bounds=([.7, -1, 0, -.075], [1.3, 1, .55, .075]),
                            loss="soft_l1", f_scale=.012)
        depth, shift = fit.x[2:]
        if fit.success and depth > .035 and abs(shift) < .074:
            yield dict(clock=name, cycle=cycle, epoch_bjd=mid+shift+2457000,
                       nominal_btjd=mid, residual_min=shift*1440,
                       fit_depth=depth, rms=np.std(fit.fun), cadences=len(x))


def measure():
    events = []
    for sector in SECTORS:
        data = pd.read_csv(HERE / "data" / f"s{sector}.csv.gz")
        time = data.time.to_numpy()
        cadence = np.median(np.diff(time))
        for aperture in ["box", "central"]:
            counts = data[aperture]
            baseline = pd.Series(counts.values, index=pd.to_timedelta(time, unit="D")).rolling(
                "0.9D", center=True, min_periods=4).median()
            flux = counts.to_numpy() / baseline.to_numpy()
            valid = np.isfinite(flux) & (baseline.to_numpy() > 0)
            for half_window in [.13, .10]:
                for name in CLOCKS:
                    for row in one_clock(time[valid], flux[valid], name, half_window, cadence):
                        events.append(dict(sector=sector, aperture=aperture,
                                           half_window_d=half_window, **row))
    return pd.DataFrame(events)


def contrast(t, y):
    early_slope = (y[1] - y[0]) / (t[1] - t[0])
    return (y[3] - y[2]) - early_slope * (t[3] - t[2])


def summarize(events):
    rng = np.random.default_rng(43472931)
    out = []
    visits = []
    for (aperture, window, clock), group in events.groupby(["aperture", "half_window_d", "clock"]):
        data = [group[group.sector == s] for s in SECTORS]
        t = np.array([v.nominal_btjd.median() for v in data])
        y = np.array([v.residual_min.median() for v in data])
        # Draw contiguous blocks of three qualifying eclipse timings. This
        # assesses sampling variability within visits; the aperture/window
        # checks address some reduction systematics separately.
        boot = []
        for v in data:
            vals = v.sort_values("cycle").residual_min.to_numpy()
            n = len(vals)
            starts = rng.integers(0, n, (4000, int(np.ceil(n / 3))))
            idx = ((starts[:, :, None] + np.arange(3)) % n).reshape(4000, -1)[:, :n]
            boot.append(np.median(vals[idx], axis=1))
        ci = np.percentile(contrast(t, np.array(boot)), [2.5, 97.5])
        out.append(dict(aperture=aperture, half_window_d=window, clock=clock,
                        change_min=contrast(t, y), block_low=ci[0], block_high=ci[1]))
        for s, v, tm, dt in zip(SECTORS, data, t, y):
            visits.append(dict(aperture=aperture, half_window_d=window, clock=clock,
                               sector=s, qualifying=len(v), median_btjd=tm, residual_min=dt))
    return pd.DataFrame(visits), pd.DataFrame(out)


if __name__ == "__main__":
    events = measure()
    visits, changes = summarize(events)
    events.to_csv(HERE / "timings.csv.gz", index=False, float_format="%.10f",
                  compression={"method": "gzip", "mtime": 0})
    visits.to_csv(HERE / "visits.csv", index=False, float_format="%.7f")
    changes.to_csv(HERE / "changes.csv", index=False, float_format="%.7f")
