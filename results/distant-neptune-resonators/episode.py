"""Check finite-time trapping around the present, including phase turns and circulation."""

from concurrent.futures import ProcessPoolExecutor
import csv
import numpy as np
from scipy.signal import find_peaks
from clones import setup
from dynamics import WORK, integrate, angle

TARGET = {"spkid": "20445473", "p": 23, "q": 2}


def compute(spec):
    mode, sign, n, years, step = spec
    sim = setup(TARGET, n, sign * step)
    sim.integrator = mode
    times, data = integrate(sim, sign * years, 100)
    traces = np.array(
        [np.rad2deg(np.unwrap(angle(data[:, k], 23, 2))) for k in range(n + 1)]
    )
    return mode, sign, times, traces


if __name__ == "__main__":
    specs = [("whfast", i, 16, 5000000, 0.5) for i in (-1, 1)] + [
        ("ias15", 1, 0, 100000, 0.5)
    ]
    with ProcessPoolExecutor(3) as pool:
        results = list(pool.map(compute, specs))
    _, _, tm, back = results[0]
    _, _, tp, fwd = results[1]
    _, _, ti, ias = results[2]
    with (WORK / "episode.csv").open("w") as handle:
        wr = csv.writer(handle)
        wr.writerow(
            [
                "clone",
                "range_plusminus_0p5Myr_deg",
                "range_plusminus_1Myr_deg",
                "past_first_360yr",
                "future_first_360yr",
                "past_median_peak_period_yr",
                "future_median_peak_period_yr",
            ]
        )
        for i, (b, f) in enumerate(zip(back, fwd)):

            def duration(t, trace):
                crossed = np.flatnonzero(
                    np.maximum.accumulate(trace) - np.minimum.accumulate(trace) >= 360
                )
                return abs(int(t[crossed[0]])) if len(crossed) else ""

            def period(t, trace):
                series = trace[np.abs(t) <= 400000]
                peaks, _ = find_peaks(series, prominence=40, distance=80)
                return int(np.median(np.diff(peaks)) * 100) if len(peaks) > 1 else ""

            both = np.r_[b[1:10001][::-1], f[:10001]]
            wr.writerow(
                [
                    i,
                    round(np.ptp(np.r_[b[:5001], f[:5001]]), 3),
                    round(np.ptp(both), 3),
                    duration(tm, b),
                    duration(tp, f),
                    period(tm, b),
                    period(tp, f),
                ]
            )
    with (WORK / "nominal-track.csv").open("w") as handle:
        wr = csv.writer(handle)
        wr.writerow(["offset_years", "whfast_phi_deg", "ias15_phi_deg"])
        t = np.r_[tm[1:][::-1], tp]
        phi = np.r_[back[0, 1:][::-1], fwd[0]]
        for yr, angle_deg in zip(t[::10], phi[::10]):
            if abs(yr) <= 1100000:
                wr.writerow(
                    [
                        int(yr),
                        round(angle_deg % 360, 4),
                        round(ias[0, int(yr / 100)] % 360, 4)
                        if 0 <= yr <= ti[-1]
                        else "",
                    ]
                )
    dif = (ias[0] - fwd[0, : len(ti)] + 180) % 360 - 180
    print("max 100kyr IAS15 vs WHFast phase difference [deg]", np.max(abs(dif)))
