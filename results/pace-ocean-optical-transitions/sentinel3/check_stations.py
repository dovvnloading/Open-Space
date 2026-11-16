"""Independent OLCI spot, spatial and seasonal checks for additional same-station candidates."""

import datetime as dt
import itertools
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from olci import distance, save
from olci_zarr import RECEIPTS, Grid

# Dates and index choices selected from the independent-of-OLCI PACE/field screen.
CANDIDATES = [
    ("Casablanca_Platform", 40.717, 1.358, "2024-08-01", "2024-08-02", (560, 620, 665)),
    ("AAOT", 45.3139, 12.5083, "2025-11-06", "2025-11-07", (560, 620, 665)),
    (
        "Section-7_Platform",
        44.5458,
        29.4466,
        "2026-03-09",
        "2026-03-12",
        (443, 490, 560),
    ),
    (
        "Irbe_Lighthouse",
        57.75092,
        21.72297,
        "2024-06-20",
        "2024-06-24",
        (443, 490, 560),
    ),
    ("USC_SEAPRISM", 33.56371, -118.11782, "2025-08-11", "2025-08-16", (443, 490, 560)),
]


def shape(rrs, bands):
    a, b, c = bands
    first, middle, last = rrs
    with np.errstate(divide="ignore", invalid="ignore"):
        residual = np.log(middle) - (
            (c - b) * np.log(first) + (b - a) * np.log(last)
        ) / (c - a)
    return np.where((first > 0) & (middle > 0) & (last > 0), residual, np.nan)


def run():
    event_grid, history = Grid("timeChunked"), Grid("geoChunked")
    station_rows, ref_rows = [], []
    for station, lat, lon, start, end, bands in CANDIDATES:
        wave = [f"RRS{w}" for w in bands]
        ev, finish = map(dt.date.fromisoformat, (start, end))
        lag = finish - ev
        y, x = event_grid.bounds(lat - 1.1, lat + 1.1, lon - 2, lon + 2)
        ys, xs = event_grid.lat[y], event_grid.lon[x]
        near = distance(ys[:, None], xs[None, :], lat, lon)
        i, j = np.abs(ys - lat).argmin(), np.abs(xs - lon).argmin()

        def pull(k, y=y, x=x):
            d, w = k
            return event_grid.extract(w, d, y, x)

        with ThreadPoolExecutor(max_workers=9) as pool:
            arrays = list(pool.map(pull, itertools.product((ev, finish), wave)))
        delt = shape(arrays[3:], bands) - shape(arrays[:3], bands)
        ix, jx = np.abs(history.lat - lat).argmin(), np.abs(history.lon - lon).argmin()
        iy, iz = np.arange(ix - 1, ix + 2), np.arange(jx - 1, jx + 2)
        by, bz = (
            np.unique(iy // history.meta[wave[0] + "/.zarray"]["chunks"][1]),
            np.unique(iz // history.meta[wave[0] + "/.zarray"]["chunks"][2]),
        )
        nt, ny, nx = history.meta[wave[0] + "/.zarray"]["chunks"]
        dates = sorted(history.dates, key=history.dates.get)
        times = [
            (d, d + lag)
            for d in dates
            if d.year != ev.year
            and abs((d - ev.replace(year=d.year)).days) <= 15
            and d + lag in history.dates
        ]
        blocks = sorted({history.dates[d] // nt for pair in times for d in pair})
        data = np.full((3, len(dates), 3, 3), np.nan)
        jobs = list(itertools.product(wave, blocks, by, bz))

        def block(k, iy=iy, iz=iz, ny=ny, nx=nx):
            w, t, k1, k2 = k
            return history.chunk(w, [t, k1, k2])[:, iy[iy // ny == k1] % ny, :][
                :, :, iz[iz // nx == k2] % nx
            ]

        with ThreadPoolExecutor(max_workers=10) as pool:
            for (w, t, k1, k2), a in zip(jobs, pool.map(block, jobs)):
                index = np.arange(t * nt, min((t + 1) * nt, len(dates)))
                dest = np.ix_(
                    index,
                    np.flatnonzero(iy // ny == k1),
                    np.flatnonzero(iz // nx == k2),
                )
                data[wave.index(w)][dest] = a[: len(index)]
        curve = shape(data, bands)
        reference = []
        for a, b in times:
            z = curve[history.dates[b]] - curve[history.dates[a]]
            valid = np.isfinite(z)
            if valid.sum() >= 5:
                med = np.median(z[valid])
                reference.append(med)
                ref_rows.append([station, start, end, a, b, valid.sum(), med, z[1, 1]])
        reference = np.array(reference)
        nhood = delt[i - 1 : i + 2, j - 1 : j + 2]
        finite = nhood[np.isfinite(nhood)]
        med = np.median(finite) if finite.size else np.nan
        n_extreme = np.sum(reference <= med) if med < 0 else np.sum(reference >= med)
        for rad in (10, 25, 50, 100):
            area = near <= rad
            valid = delt[area & np.isfinite(delt)]
            station_rows.append(
                [
                    station,
                    lat,
                    lon,
                    start,
                    end,
                    *bands,
                    ys[i],
                    xs[j],
                    delt[i, j],
                    finite.size,
                    med,
                    len(reference),
                    n_extreme if finite.size >= 5 else "",
                    rad,
                    area.sum(),
                    len(valid),
                    np.median(valid) if valid.size else "",
                    np.mean(valid < 0) if valid.size else "",
                ]
            )
        print(
            station,
            "center",
            delt[i, j],
            "paired",
            finite.size,
            "median",
            med,
            "history",
            len(reference),
            "as extreme",
            n_extreme,
            "50 km",
            station_rows[-2][-3:],
            flush=True,
        )
    save(
        "other_station_checks.csv",
        "station,lat,lon,start,end,first_nm,center_nm,last_nm,olci_lat,olci_lon,center_change,paired_3x3,median_change_3x3,historical_comparisons,historical_same_direction_as_extreme,radius_km,grid_cells,paired_radius_cells,radius_median_change,radius_fraction_negative",
        station_rows,
    )
    save(
        "other_station_reference.csv",
        "station,event_start,event_end,reference_start,reference_end,paired_cells,median_change_3x3,center_change",
        ref_rows,
    )
    save(
        "other_station_receipts.csv",
        "url,bytes,sha256",
        [(u, *v) for u, v in sorted(RECEIPTS.items())],
    )


if __name__ == "__main__":
    run()
