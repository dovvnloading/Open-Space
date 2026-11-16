"""Compare OLCI spectral shape on the event footprint against preceding days and historical analogs."""

import csv
import datetime as dt
import itertools
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from check_stations import shape
from olci import save
from olci_zarr import RECEIPTS, ROOT, Grid

START = dt.date(2025, 11, 6)
END = dt.date(2025, 11, 7)
BANDS = [443, 490, 510, 560, 620, 665, 674, 681, 709]


def main():
    with (ROOT / "aaot_paired_spectra.csv").open() as f:
        target = [
            r
            for r in csv.DictReader(f)
            if float(r["distance_km"]) <= 50 and float(r["depth_m"]) > 0
        ]
    with (ROOT / "aaot_regional_seasonal.csv").open() as f:
        analogs = [
            (
                dt.date.fromisoformat(r["reference_start"]),
                dt.date.fromisoformat(r["reference_end"]),
            )
            for r in csv.DictReader(f)
            if float(r["historical_median"]) <= float(r["matched_event_median"])
        ]
    pairs = [
        (START - dt.timedelta(days=n), END - dt.timedelta(days=n)) for n in (3, 2, 1, 0)
    ] + analogs
    pairs += [(START - dt.timedelta(days=1), END)]
    days = sorted(set(itertools.chain.from_iterable(pairs)))
    grid = Grid("timeChunked")
    y, x = grid.bounds(44.7, 46.0, 11.5, 13.5)
    ys, xs = grid.lat[y], grid.lon[x]
    rows = np.array([np.abs(ys - float(r["lat"])).argmin() for r in target])
    cols = np.array([np.abs(xs - float(r["lon"])).argmin() for r in target])
    area = np.array(
        [
            6371**2
            * np.radians(1 / 24)
            * (
                np.sin(np.radians(float(r["lat"]) + 1 / 48))
                - np.sin(np.radians(float(r["lat"]) - 1 / 48))
            )
            for r in target
        ]
    )
    keys = list(itertools.product(days, BANDS))

    def load(k):
        day, w = k
        return grid.extract(f"RRS{w}", day, y, x)[rows, cols]

    with ThreadPoolExecutor(max_workers=12) as pool:
        data = dict(zip(keys, pool.map(load, keys)))
    out = []
    for a, b in pairs:
        valid = np.ones(len(rows), dtype=bool)
        for day, w in itertools.product((a, b), BANDS):
            valid &= np.isfinite(data[day, w]) & (data[day, w] > 0)
        z = []
        for day in (a, b):
            # The linear continuum at 681 uses the measured 665 and 709 nm as anchors.
            z.append(data[day, 681] - (28 * data[day, 665] + 16 * data[day, 709]) / 44)
        if valid.any():
            vals = [
                np.median(np.log(data[b, w][valid] / data[a, w][valid])) for w in BANDS
            ]
            q = []
            for ws in [(443, 490, 560), (560, 620, 665), (665, 681, 709)]:
                t0 = shape([data[a, w] for w in ws], ws)
                t1 = shape([data[b, w] for w in ws], ws)
                q.append(np.median((t1 - t0)[valid]))
            fl1, fl2 = np.median(z[0][valid]), np.median(z[1][valid])
            difference = np.median((z[1] - z[0])[valid])
        else:
            vals = [""] * len(BANDS)
            q = [""] * 3
            fl1 = fl2 = difference = ""
        out.append(
            [
                a,
                b,
                len(target),
                valid.sum(),
                area[valid].sum(),
                *vals,
                *q,
                fl1,
                fl2,
                difference,
            ]
        )
    save(
        "aaot_spectral_context.csv",
        "start,end,event_footprint_cells,paired_all_bands,paired_area_km2,"
        + ",".join(f"median_log_change_{w}" for w in BANDS)
        + ",median_delta_C490,median_delta_C620,median_delta_C681,median_start_LH681_sr-1,median_end_LH681_sr-1,median_delta_LH681_sr-1",
        out,
    )
    save(
        "aaot_spectral_receipts.csv",
        "url,bytes,sha256",
        [(u, *v) for u, v in sorted(RECEIPTS.items())],
    )


if __name__ == "__main__":
    main()
