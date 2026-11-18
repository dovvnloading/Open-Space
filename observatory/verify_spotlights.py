"""Check the numbers printed on the observing deck against a single revision's tables."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
import csv
import io
import json
import math
import re
from statistics import correlation


def rows(read_blob, path):
    return list(csv.DictReader(io.StringIO(read_blob(path).decode("utf-8"))))


def separation(ra1, dec1, ra2, dec2):
    """Great-circle distance in arcseconds for two catalog positions in degrees."""
    a, b = map(math.radians, (dec1, dec2))
    h = math.sin((a - b) / 2) ** 2
    h += math.cos(a) * math.cos(b) * math.sin(math.radians(ra2 - ra1) / 2) ** 2
    return math.degrees(2 * math.asin(math.sqrt(h))) * 3600


def mars(read_blob):
    folder = "results/maven-2025-sep-events/"
    table, survey = folder + "daily.csv", folder + "summary.json"
    days = [row for row in rows(read_blob, table) if row["date"].startswith("2025-08-")]
    summary = json.loads(read_blob(survey))
    event, = [e for e in summary["ranked_events"] if e["start_utc"].startswith("2025-08-22")]
    result = {}
    for species in ("ion", "electron"):
        field = f"max_hourly_median_{species}"
        peak = max(float(day[field]) for day in days)
        if peak != event[f"peak_hourly_median_{species}"]:
            raise ValueError(f"MAVEN {species} peak disagrees between the daily and event tables")
        reference = summary["survey_hourly_reference"][f"median_{species}_flux_at_100plus_sample_hours"]
        value = f"≈{round(peak / reference, -1):,.0f}×"
        method = (f"Divide the August peak hourly median ({peak:,.0f}) by the survey median "
                  f"({reference:.6g}) for hours with at least 100 samples. Both fluxes are in (cm² s sr)⁻¹.")
        result[species] = value, method, [table, survey]
    result["orbits"] = (str(event["distinct_joint_orbits"]),
                        "Distinct spacecraft orbits with joint ion and electron samples in the selected hours.",
                        [survey])
    table_view = {
        "caption": "Every observed UTC date in August 2025. Ions and electrons are the highest hourly "
                   "median fluxes for each date, in (cm² s sr)⁻¹. Values retain the daily table’s precision.",
        "columns": ["UTC date", "Observed hours", "Selected hours", "Ions", "Electrons"],
        "rows": [[r[k] for k in ("date", "observed_hours", "event_hours",
                                 "max_hourly_median_ion", "max_hourly_median_electron")] for r in days],
        "sourcePath": table,
    }
    return result, table_view


def sun(read_blob):
    folder = "results/2021-july-13-farside-eruption/inputs/"
    positions = folder + "stix_extract.csv"
    summary = folder + "noaa_srs_20210722.txt"
    flare, = [r for r in rows(read_blob, positions) if r["flare_id"] == "2107131835"]
    measured = math.degrees(math.atan2(-float(flare["hpc_x_solo"]), float(flare["hpc_y_solo"]))) % 360
    result = {"angle": (f"{measured:.2f}°", "Position angle east from solar north, calculated from "
                        "the published STIX image centroid for flare 2107131835.", [positions])}
    longitude = sum(float(flare[k]) * sign for k, sign in
                    (("hgc_lon", 1), ("hgs_lon", -1), ("solo_position_lon", 1)))
    observer_lat = math.radians(float(flare["solo_position_lat"]))
    regions = re.findall(r"^(\d{4})\s+([NS])(\d+)[EW]\d+\s+(\d+)\b",
                         read_blob(summary).decode(), re.MULTILINE)
    residuals = []
    for region, hemisphere, latitude, carrington in regions:
        lat = math.radians(int(latitude) * (-1 if hemisphere == "S" else 1))
        delta = math.radians(int(carrington) - longitude)
        x = math.cos(lat) * math.sin(delta)
        y = math.sin(lat) * math.cos(observer_lat) - math.cos(lat) * math.cos(delta) * math.sin(observer_lat)
        predicted = math.degrees(math.atan2(-x, y))
        residuals.append((abs((measured - predicted + 180) % 360 - 180), region))
    if [region for _, region in sorted(residuals)[:2]] != ["2849", "2847"]:
        raise ValueError("STIX candidate ordering differs from the two regions discussed in the dossier")
    for residual, region in sorted(residuals)[:2]:
        result[region] = (f"{residual:.2f}°", f"Project NOAA region 1{region} from the 22 July summary "
                          "into Solar Orbiter’s 13 July view and compare its direction with the X-ray centroid.",
                          [positions, summary])
    rates_path = folder + "stix_rates.csv"
    rates = [r for r in rows(read_blob, rates_path) if r["product"] == "SCI" and r["band_keV"] == "15-25"
             and "2021-07-13T18:15" <= r["spacecraft_utc"] < "2021-07-13T18:45"]
    if not rates:
        raise ValueError("No STIX science accumulations in the displayed hard X-ray interval")
    table_view = {
        "caption": "All 15–25 keV science-product accumulations from 18:15 up to 18:45 on 13 July 2021. "
                   "The times are measured at the spacecraft. Count rates are rounded to three decimals.",
        "columns": ["Spacecraft UTC", "Count rate (counts/s)"],
        "rows": [[r["spacecraft_utc"].split("T")[1], f"{float(r['count_rate_s']):.3f}"] for r in rates],
        "sourcePath": rates_path,
    }
    return result, table_view


def smooth(samples, grid, width):
    """The released lens scan uses a two-width window and an error floor of 0.035 mag."""
    dates = [sample[0] for sample in samples]
    smoothed = []
    for day in grid:
        left = bisect_right(dates, day - 2 * width)
        right = bisect_left(dates, day + 2 * width)
        weight_sum = weighted_mag = 0.0
        for t, mag, error in samples[left:right]:
            weight = math.exp(-0.5 * ((day - t) / width) ** 2) / (error ** 2 + 0.035 ** 2)
            weight_sum += weight
            weighted_mag += weight * mag
        smoothed.append(weighted_mag / weight_sum if weight_sum else math.nan)
    return smoothed


def peak_delay(table, time_column, band):
    """Independent scalar replay of the published 5-day correlation scan, B leading A."""
    curves = {image: sorted((float(r[time_column]), float(r["magnitude"]), float(r["median_catalog_error"]))
                            for r in table if r["image"] == image and r["band"] == band)
              for image in ("A", "B")}
    dates = [t for curve in curves.values() for t, _, _ in curve]
    grid = range(5 * (math.floor(min(dates) / 5) - 1), 5 * (math.floor(max(dates) / 5) + 2), 5)
    residuals = {image: [x - trend for x, trend in zip(smooth(samples, grid, 10), smooth(samples, grid, 140))]
                 for image, samples in curves.items()}
    scores = []
    for lag in range(-600, 601, 5):
        pairs = []
        for i, leading in enumerate(residuals["B"]):
            other = i + lag // 5
            if 0 < other < len(grid) - 1 and math.isfinite(leading):
                following = residuals["A"][other]
                if math.isfinite(following):
                    pairs.append((leading, following))
        if len(pairs) > 20:
            scores.append((lag, correlation(*zip(*pairs))))
    if not scores:
        raise ValueError("Lens scan has insufficient overlapping nights")
    return max(scores, key=lambda pair: pair[1])[0]


def quasar(read_blob):
    folder = "results/j0826-lens-delay/"
    ztf_path, gaia_path, pos_path = [folder + name for name in ("ztf-nights.csv", "gaia-nights.csv", "targets.csv")]
    ztf = [r for r in rows(read_blob, ztf_path) if r["system"] == "J0826+7002"]
    gaia = rows(read_blob, gaia_path)
    pos = {r["image"]: r for r in rows(read_blob, pos_path) if r["system"] == "J0826+7002"}
    arcsec = separation(*(float(pos[image][key]) for image in ("A", "B") for key in ("ra_deg", "dec_deg")))
    measured = {
        "ztf": (f"{peak_delay(ztf, 'mjd', 'zg')} days", "Peak of the detrended ZTF g-band correlation "
                "on a five-day grid, with image B leading A.", [ztf_path]),
        "gaia": (f"{peak_delay(gaia, 'bjd_tcb_minus_2400000p5', 'g')} days",
                 "Same five-day search on independent Gaia G-band nights, with image B leading A.", [gaia_path]),
        "separation": (f"{arcsec:.2f}″", "Great-circle distance between the listed A and B image coordinates.",
                       [pos_path]),
    }
    time = "bjd_tcb_minus_2400000p5"
    g_nights = sorted((r for r in gaia if r["band"] == "g"), key=lambda r: (float(r[time]), r["image"]))
    table_view = {
        "caption": "All accepted Gaia G-band nights for images A and B. Each row retains its native epoch "
                   "and measured magnitude before the horizontal and vertical offsets used in the figure.",
        "columns": ["Image", "BJD TCB − 2,400,000.5", "G (mag)", "Median catalog error (mag)"],
        "rows": [[r["image"], r[time], r["magnitude"], r["median_catalog_error"]] for r in g_nights],
        "sourcePath": gaia_path,
    }
    return measured, table_view


def dwarf(read_blob):
    path = "results/wise-j0458-radio/epochs.csv"
    source_rows = rows(read_blob, path)
    epochs = {int(r["epoch"]): r for r in source_rows}
    if len(source_rows) != 4 or set(epochs) != {1, 2, 3, 4}:
        raise ValueError("Expected one row per VLASS visit in the four-epoch study")

    def displacement(prefix):
        return separation(*(float(epochs[i][prefix + key]) for i in (1, 4) for key in ("ra", "dec")))

    flux = float(epochs[2]["forced_mjy"])
    noise = float(epochs[2]["control_mad_mjy"])
    measured = {
        "radio": (f"{displacement('centroid_'):.2f}″", "Great-circle distance between the freely fitted radio "
                  "centroids in October 2017 and September 2025.", [path]),
        "infrared": (f"{displacement(''):.2f}″", "Displacement of the catalog proper-motion prediction "
                     "between the same two radio visits.", [path]),
        "nondetection": (f"{flux:.2f}".replace("-", "−") + " mJy/beam",
                         f"Beam fit at the predicted position in August 2020. Nearby control fits "
                         f"give a scatter of {noise:.2f} mJy/beam.", [path]),
    }
    table_rows = []
    for _, visit in sorted(epochs.items()):
        centroid = "No centroid recovered"
        if visit["centroid_ra"] and visit["centroid_dec"]:
            centroid = f"{float(visit['centroid_ra']):.6f}, {float(visit['centroid_dec']):+.6f}"
        table_rows.append([visit["date"].split("T")[0], f"{float(visit['forced_mjy']):+.3f}",
                           f"{float(visit['control_mad_mjy']):.3f}", centroid])
    table_view = {
        "caption": "All four VLASS visits. Flux is fitted at the predicted infrared position. Scatter comes "
                   "from fits in nearby sky. The independent centroid fit is printed to six decimal places.",
        "columns": ["UTC date", "Flux (mJy/beam)", "Scatter (mJy/beam)", "Centroid RA, Dec (°)"],
        "rows": table_rows,
        "sourcePath": path,
    }
    return measured, table_view


def verify_spotlights(items, read_blob):
    """Fail on missing, duplicate or drifting readouts. Attach the paths behind each method."""
    methods = {"mars": mars, "sun": sun, "quasar": quasar, "dwarf": dwarf}
    ids = [item["id"] for item in items]
    if len(ids) != len(set(ids)) or set(ids) != set(methods):
        raise ValueError("Expected one dossier for each audited target")
    for item in items:
        expected, item["table"] = methods[item["id"]](read_blob)
        metrics = item["measurements"]
        if len(metrics) != len(expected) or {m.get("key") for m in metrics} != set(expected):
            raise ValueError(f"Measurement keys differ from the audit for {item['id']}")
        for metric in metrics:
            value, method, paths = expected[metric["key"]]
            if value != metric["value"]:
                raise ValueError(f"{item['id']}.{metric['key']}: displayed {metric['value']!r}, "
                                 f"source tables give {value!r}")
            metric.update(method=method, sourcePaths=paths)
