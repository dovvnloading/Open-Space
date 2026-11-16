"""Verify endpoint reflectance, regional extent and historical context for the optical screen."""

import csv
import datetime as dt
import json
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode

import numpy as np
from olci_zarr import RECEIPTS, ROOT, Grid, read

WAVES = ["RRS443", "RRS490", "RRS560", "RRS665"]
USC = (33.56371, -118.11782)
EVENTS = [("2024-03-11", "2024-03-15"), ("2024-04-10", "2024-04-15")]
GULF = [
    (27.6981, -83.2478, "2024-05-22"),
    (27.5157, -83.3684, "2024-05-22"),
    (24.7922, -84.0410, "2024-05-24"),
    (24.7917, -83.8267, "2024-05-24"),
    (24.7866, -83.7095, "2024-05-24"),
    (27.9120, -83.2203, "2024-05-22"),
    (24.7956, -83.9243, "2024-05-24"),
]


def save(name, header, rows):
    with (ROOT / name).open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header.split(","))
        writer.writerows(rows)


def curvature(rrs):
    blue, cyan, green = rrs[:3]
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.log(cyan) - 70 / 117 * np.log(blue) - 47 / 117 * np.log(green)
    return np.where((blue > 0) & (cyan > 0) & (green > 0), c, np.nan)


def distance(lat, lon, target_lat, target_lon):
    a = np.sin(np.radians(lat - target_lat) / 2) ** 2
    a = (
        a
        + np.cos(np.radians(lat))
        * np.cos(np.radians(target_lat))
        * np.sin(np.radians(lon - target_lon) / 2) ** 2
    )
    return 6371 * 2 * np.arcsin(np.sqrt(a))


def extract_events():
    grid = Grid("timeChunked")
    dates = [day for pair in EVENTS for day in pair]
    regions = {day: grid.bounds(32, 35, -120, -116) for day in dates}
    regions.update({day: grid.bounds(24.5, 28.2, -84.5, -82.9) for _, _, day in GULF})
    jobs = [(day, wave) for day in regions for wave in WAVES]

    def load(job):
        day, wave = job
        return grid.extract(wave, dt.date.fromisoformat(day), *regions[day])

    with ThreadPoolExecutor(max_workers=10) as pool:
        bands = dict(zip(jobs, pool.map(load, jobs)))
    spectra = {day: np.stack([bands[day, wave] for wave in WAVES]) for day in regions}
    curves = {day: curvature(rrs) for day, rrs in spectra.items()}
    stations = []
    for lat, lon, day in [(*USC, day) for day in dates] + GULF:
        y, x = regions[day]
        ys, xs = grid.lat[y], grid.lon[x]
        i, j = np.abs(ys - lat).argmin(), np.abs(xs - lon).argmin()
        win = curves[day][i - 1 : i + 2, j - 1 : j + 2]
        valid = win[np.isfinite(win)]
        stats = (
            [np.median(valid), valid.min(), valid.max()] if valid.size else ["", "", ""]
        )
        stations.append(
            [
                day,
                lat,
                lon,
                ys[i],
                xs[j],
                distance(lat, lon, ys[i], xs[j]),
                *spectra[day][:, i, j],
                curves[day][i, j],
                valid.size,
                *stats,
            ]
        )
    save(
        "olci_stations.csv",
        "date,station_lat,station_lon,olci_lat,olci_lon,distance_km,RRS443,RRS490,RRS560,RRS665,C,valid_3x3,C_median_3x3,C_min_3x3,C_max_3x3",
        stations,
    )
    spatial = []
    paired = []
    for start, end in EVENTS:
        y, x = regions[start]
        ys, xs = grid.lat[y], grid.lon[x]
        dist = distance(ys[:, None], xs[None, :], *USC)
        change = curves[end] - curves[start]
        lat_bound = np.sin(np.radians(ys.astype(float) + 1 / 48)) - np.sin(
            np.radians(ys.astype(float) - 1 / 48)
        )
        area = np.broadcast_to(
            6371**2 * np.radians(1 / 24) * lat_bound[:, None], dist.shape
        )
        for radius in (10, 25, 50, 100, 125):
            within = dist <= radius
            valid = within & np.isfinite(change)
            delta = change[valid]
            qs = np.percentile(delta, [50, 10, 90]) if delta.size else [""] * 3
            fractions = [
                a.mean() if delta.size else ""
                for a in (delta < 0, delta > 0, delta <= -0.25, delta >= 0.25)
            ]
            spatial.append(
                [
                    start,
                    end,
                    radius,
                    within.sum(),
                    delta.size,
                    area[valid].sum(),
                    *qs,
                    *fractions,
                ]
            )
        for i, j in zip(*np.where((dist <= 125) & np.isfinite(change))):
            paired.append(
                [
                    start,
                    end,
                    ys[i],
                    xs[j],
                    dist[i, j],
                    area[i, j],
                    curves[start][i, j],
                    curves[end][i, j],
                    change[i, j],
                ]
            )
    save(
        "olci_spatial.csv",
        "date_before,date_after,radius_km,grid_cells,paired_cells,area_km2,delta_C_median,delta_C_p10,delta_C_p90,fraction_declining,fraction_increasing,fraction_change_le_minus_0.25,fraction_change_ge_0.25",
        spatial,
    )
    save(
        "olci_paired_cells.csv",
        "date_before,date_after,lat,lon,distance_km,cell_area_km2,C_before,C_after,delta_C",
        paired,
    )


def seasonal_reference():
    grid = Grid("geoChunked")
    dates = sorted(grid.dates, key=grid.dates.get)
    i = np.abs(grid.lat - USC[0]).argmin()
    j = np.abs(grid.lon - USC[1]).argmin()
    y, x = np.arange(i - 1, i + 2), np.arange(j - 1, j + 2)
    nt, ny, nx = grid.meta["RRS443/.zarray"]["chunks"]
    jobs = [
        (wave, t, a, b)
        for wave in WAVES[:3]
        for t in range((len(dates) - 1) // nt + 1)
        for a in np.unique(y // ny)
        for b in np.unique(x // nx)
    ]
    series = np.full((3, len(dates), 3, 3), np.nan)

    def load(job):
        wave, t, a, b = job
        return grid.chunk(wave, [t, a, b])[:, y[y // ny == a] % ny, :][
            :, :, x[x // nx == b] % nx
        ]

    with ThreadPoolExecutor(max_workers=10) as pool:
        for (wave, t, a, b), block in zip(jobs, pool.map(load, jobs)):
            ix = np.arange(t * nt, min((t + 1) * nt, len(dates)))
            selection = np.ix_(
                ix, np.flatnonzero(y // ny == a), np.flatnonzero(x // nx == b)
            )
            series[WAVES.index(wave)][selection] = block[: len(ix)]
    curves = curvature(series)
    references, summary = [], []
    for start, end in EVENTS:
        ev, last = map(dt.date.fromisoformat, (start, end))
        lag = last - ev
        target = curves[grid.dates[last]] - curves[grid.dates[ev]]
        paired = target[np.isfinite(target)]
        rows = []
        for d in dates:
            if (
                d.year == ev.year
                or abs((d - ev.replace(year=d.year)).days) > 15
                or d + lag not in grid.dates
            ):
                continue
            delta = curves[grid.dates[d + lag]] - curves[grid.dates[d]]
            valid = delta[np.isfinite(delta)]
            if valid.size >= 5:
                rows.append(
                    [
                        start,
                        d,
                        d + lag,
                        lag.days,
                        valid.size,
                        np.median(valid),
                        delta[1, 1],
                    ]
                )
        references.extend(rows)
        local = np.array([r[-2] for r in rows])
        center = np.array([r[-1] for r in rows])
        center = center[np.isfinite(center)]
        summary.append(
            [
                start,
                end,
                paired.size,
                np.median(paired),
                target[1, 1],
                len(local),
                sum(local <= np.median(paired)),
                len(center),
                sum(center <= target[1, 1]),
            ]
        )
    save(
        "olci_seasonal_reference.csv",
        "event_date_before,reference_date_before,reference_date_after,lag_days,paired_pixels,median_same_pixel_change_3x3,center_pixel_change",
        references,
    )
    save(
        "olci_seasonal_summary.csv",
        "date_before,date_after,paired_event_pixels,event_median_change_3x3,event_center_change,neighborhood_comparisons,neighborhood_as_negative,center_comparisons,center_as_negative",
        summary,
    )


def runoff():
    parameters = {
        "f": "json",
        "monitoring_location_id": "USGS-11078000",
        "parameter_code": "00060",
        "statistic_id": "00003",
        "datetime": "2024-03-01/2024-04-22",
        "limit": 1000,
    }
    url = (
        "https://api.waterdata.usgs.gov/ogcapi/v0/collections/daily/items?"
        + urlencode(parameters, safe="/")
    )
    features = json.loads(read(url))["features"]
    fields = [
        "time",
        "monitoring_location_id",
        "value",
        "unit_of_measure",
        "approval_status",
        "qualifier",
    ]
    rows = sorted([[f["properties"][k] for k in fields] for f in features])
    if len(rows) != 53 or len({r[0] for r in rows}) != 53:
        raise ValueError("Unexpected number of station days")
    save("usgs_discharge.csv", ",".join(fields), rows)


if __name__ == "__main__":
    extract_events()
    seasonal_reference()
    runoff()
    save(
        "source_receipts.csv",
        "url,bytes,sha256",
        [(u, *v) for u, v in sorted(RECEIPTS.items())],
    )
