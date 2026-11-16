"""Evaluate spatial extent, bottom depth, red brightness and prior November variability at AAOT."""

import datetime as dt
import io
import itertools
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import quote

import numpy as np
from check_stations import shape
from olci import distance, save
from olci_zarr import RECEIPTS, Grid, read
from scipy.io import netcdf_file

LAT, LON = 45.3139, 12.5083
START, END = dt.date(2025, 11, 6), dt.date(2025, 11, 7)
DEPTH_URL = "https://upwell.pfeg.noaa.gov/erddap/griddap/etopo180.nc?" + quote(
    "altitude[(44):1:(46.6)][(10.5):1:(14.5)]", safe=",:()"
)


def main():
    grid = Grid("timeChunked")
    history = Grid("geoChunked")
    y, x = grid.bounds(44.3, 46.4, 11.1, 13.9)
    ys, xs = grid.lat[y], grid.lon[x]
    d = distance(ys[:, None], xs[None, :], LAT, LON)
    with netcdf_file(io.BytesIO(read(DEPTH_URL)), mmap=False) as ds:
        yy, xx = [ds.variables[k].data.copy() for k in ("latitude", "longitude")]
        iy = np.abs(yy[:, None] - ys).argmin(0)
        ix = np.abs(xx[:, None] - xs).argmin(0)
        depth = -ds.variables["altitude"].data[np.ix_(iy, ix)].astype(float)
    wavelengths = (510, 560, 620, 665, 674, 681, 709)
    keys = list(itertools.product((START, END), wavelengths))

    def load(k):
        day, w = k
        return grid.extract(f"RRS{w}", day, y, x)

    with ThreadPoolExecutor(max_workers=12) as pool:
        data = dict(zip(keys, pool.map(load, keys)))
    waves = (560, 620, 665)
    curves = []
    for day in (START, END):
        curves.append(shape([data[day, w] for w in waves], waves))
    delta = curves[1] - curves[0]
    positive = np.isfinite(delta)
    results = []
    r0, r1 = data[START, 665], data[END, 665]
    for radius in (25, 50, 100):
        near = d <= radius
        region = near & (depth > 0) & positive
        red = r0[region]
        q = np.percentile(red, [25, 50, 75])
        ratio = r1 / r0
        masks = {
            "all_valid": np.ones_like(positive),
            "depth_lt25": (depth < 25),
            "depth_ge25": (depth >= 25),
            "depth_ge50": (depth >= 50),
            "depth_ge100": (depth >= 100),
            "start_red_lowest_quartile": r0 <= q[0],
            "start_red_highest_quartile": r0 >= q[2],
            "red_stable_within20pct": (ratio >= 0.8) & (ratio <= 1.2),
            "endpoint_red_lt0005": (r0 < 0.0005) & (r1 < 0.0005),
        }
        for label, mask in masks.items():
            ids = region & mask
            values = delta[ids]
            p = np.percentile(values, [25, 50, 75]) if values.size else ["", "", ""]
            results.append(
                [
                    radius,
                    label,
                    np.sum(near & (depth > 0) & mask),
                    values.size,
                    *p,
                    np.mean(values < 0) if values.size else "",
                    np.median(r0[ids]) if values.size else "",
                    np.median(r1[ids]) if values.size else "",
                ]
            )
    save(
        "aaot_masks.csv",
        "radius_km,subset,mask_cells,paired_cells,delta_q25,delta_median,delta_q75,fraction_negative,median_start_rrs665,median_end_rrs665",
        results,
    )
    out = []
    for ij in zip(*np.where((d <= 100) & positive)):
        out.append(
            [
                ys[ij[0]],
                xs[ij[1]],
                d[ij],
                depth[ij],
                curves[0][ij],
                curves[1][ij],
                delta[ij],
                *[data[k][ij] for k in keys],
            ]
        )
    hdr = "lat,lon,distance_km,depth_m,C_before,C_after,delta_C," + ",".join(
        f"Rrs{b}_{day}" for day, b in keys
    )
    save("aaot_paired_spectra.csv", hdr, out)
    # Match each control to the event's actual 50-km pixel set. Compare medians on the exact same cells.
    event_mask = (d <= 50) & (depth > 0) & positive
    ey, ex = np.where(event_mask)
    iy, ix = y[ey], x[ex]
    dates = sorted(history.dates, key=history.dates.get)
    days = dt.timedelta(days=1)
    intervals = [
        (t, t + days)
        for t in dates
        if t.year != START.year
        and abs((t - START.replace(year=t.year)).days) <= 15
        and t + days in history.dates
    ]
    days_used = sorted(set(itertools.chain.from_iterable(intervals)))
    day_idx = np.array([history.dates[t] for t in days_used])
    nt, ny, nx = history.meta["RRS560/.zarray"]["chunks"]
    tiles = sorted(set(zip(iy // ny, ix // nx)))
    blocks = np.unique(day_idx // nt)
    jobs = list(itertools.product(waves, blocks, tiles))
    samples = np.full((3, len(days_used), len(iy)), np.nan)

    def block(job):
        w, t, (a, b) = job
        cols = np.flatnonzero((iy // ny == a) & (ix // nx == b))
        rows = np.flatnonzero(day_idx // nt == t)
        z = history.chunk(f"RRS{w}", [t, a, b])
        return rows, cols, z[day_idx[rows, None] % nt, iy[cols] % ny, ix[cols] % nx]

    with ThreadPoolExecutor(max_workers=12) as pool:
        for job, (rows, cols, z) in zip(jobs, pool.map(block, jobs)):
            samples[waves.index(job[0])][np.ix_(rows, cols)] = z
    hc = shape(samples, waves)
    di = {t: i for i, t in enumerate(days_used)}
    hist = []
    for first, last in intervals:
        hist_delta = hc[di[last]] - hc[di[first]]
        use = np.isfinite(hist_delta)
        n = use.sum()
        if n >= np.ceil(0.5 * len(iy)):
            match = delta[ey[use], ex[use]]
            z = hist_delta[use]
            hist.append(
                [
                    first,
                    last,
                    len(iy),
                    n,
                    np.median(match),
                    np.median(z),
                    np.percentile(z, 25),
                    np.percentile(z, 75),
                    np.mean(match < 0),
                    np.mean(z < 0),
                ]
            )
    save(
        "aaot_regional_seasonal.csv",
        "reference_start,reference_end,event_cells,paired_control_cells,matched_event_median,historical_median,historical_q25,historical_q75,matched_event_fraction_negative,historical_fraction_negative",
        hist,
    )
    print(
        "AAOT historical 50-km controls, negatives",
        sum(r[5] <= r[4] for r in hist),
        len(hist),
        "over",
        len(tiles),
        "tiles",
        flush=True,
    )
    save(
        "aaot_receipts.csv",
        "url,bytes,sha256",
        [(u, *v) for u, v in sorted(RECEIPTS.items())],
    )


if __name__ == "__main__":
    main()
