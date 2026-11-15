"""Recompute selected Venus-flyby geometry from the supplied JPL Horizons tables."""
import argparse
import csv
import gzip
import hashlib
import json
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parent
AU_KM = 149_597_870.7
VENUS_RADIUS_KM = 6051.8
CASES = {
    "2021-08-09": ("SolO", "-144", "04:05", "05:00"),
    "2023-08-21": ("PSP", "-96", "12:15", "14:55"),
    "2025-02-18": ("SolO", "-144", "20:40", "21:00"),
    "2026-12-24": ("SolO", "-144", "22:45", "23:25"),
    "2028-03-18": ("SolO", "-144", "08:05", "08:45"),
}
TARGETS = {"venus": "299", "earth": "399", "jupiter": "599", "sun": "10"}


def source_url(day, name):
    _, observer, start, stop = CASES[day]
    duration = datetime.fromisoformat(f"{day}T{stop}") - datetime.fromisoformat(f"{day}T{start}")
    step = str(int(duration.total_seconds() / 10)) if name == "venus" else "1 m"
    query = {
        "COMMAND": TARGETS[name], "CENTER": "500@" + observer, "EPHEM_TYPE": "OBSERVER",
        "START_TIME": f"{day} {start}", "STOP_TIME": f"{day} {stop}", "STEP_SIZE": step,
        "QUANTITIES": "2,13,20", "CSV_FORMAT": "YES", "CAL_FORMAT": "JD",
        "ANG_FORMAT": "DEG", "EXTRA_PREC": "YES", "TIME_DIGITS": "FRACSEC",
    }
    params = {"format": "json", **{key: f"'{value}'" for key, value in query.items()}}
    return "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(params)


def read_table(day, name, origin_jd, download):
    path = ROOT / "raw" / f"{day}_{name}.txt.gz"
    if download:
        with urllib.request.urlopen(source_url(day, name), timeout=90) as response:
            text = json.load(response)["result"]
        if "$$EOE" not in text:
            raise ValueError(f"No ephemeris returned for {day} {name}")
        path.write_bytes(gzip.compress(text.encode(), mtime=0))
    content = gzip.decompress(path.read_bytes())
    table = content.decode().split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    values = np.array([[float(row[i]) for i in (0, 3, 4, 5, 6)] for row in csv.reader(table)])
    times = (values[:, 0] - origin_jd) * 86400
    ra, dec = np.radians(values[:, 1:3]).T
    vectors = np.column_stack((np.cos(ra) * np.cos(dec), np.sin(ra) * np.cos(dec), np.sin(dec)))
    # Angular semidiameter in degrees and apparent range in AU.
    state = np.column_stack((vectors, values[:, 3] / 7200, values[:, 4]))
    receipt = [day, name, len(values), hashlib.sha256(content).hexdigest(), source_url(day, name)]
    return (times, CubicSpline(times, state, extrapolate=False), receipt)


def clearance(time, venus, source, height, source_offset):
    near, far = venus(time), source(time)
    if np.any(near[..., 4] >= far[..., 4]):
        raise ValueError("The direction being tested is closer than Venus.")
    left, right = near[..., :3], far[..., :3]
    separation = np.arctan2(np.linalg.norm(np.cross(left, right), axis=-1), np.sum(left * right, axis=-1))
    limb = np.arcsin((VENUS_RADIUS_KM + height) / (AU_KM * near[..., 4]))
    return np.degrees(separation - limb) - source_offset * far[..., 3]


def solve(day, download):
    midnight = datetime.fromisoformat(day).replace(tzinfo=timezone.utc)
    origin = midnight.timestamp() / 86400 + 2440587.5
    tables = {name: read_table(day, name, origin, download) for name in TARGETS}
    left = np.ceil(max(table[0][0] for table in tables.values()))
    right = np.floor(min(table[0][-1] for table in tables.values()))
    grid = np.arange(left, right + 1)
    rows = []
    for name in ("earth", "jupiter", "sun"):
        models = [(0, 0)]
        if name == "jupiter":
            models.extend([(0, -3), (0, 3), (100, 0), (300, 0)])
        for height, offset in models:
            def margin(time):
                return clearance(time, tables["venus"][1], tables[name][1], height, offset)

            values = margin(grid)
            crossings = np.where(np.signbit(values[1:]) != np.signbit(values[:-1]))[0]
            roots = [brentq(margin, grid[i], grid[i + 1]) for i in crossings]
            if len(roots) not in (0, 2) or values[0] < 0 or values[-1] < 0:
                raise ValueError(f"Window does not bracket all contacts: {day} {name}")
            index = int(np.argmin(values))
            minimum = values[index]
            if 0 < index < len(grid) - 1:
                minimum = minimize_scalar(margin, bounds=(grid[index - 1], grid[index + 1]), method="bounded").fun
            ingress, egress = [(midnight + timedelta(seconds=t)).isoformat(timespec="milliseconds").replace("+00:00", "Z")
                                for t in roots] if roots else ("", "")
            rows.append([day, CASES[day][0], name, height, offset, ingress, egress,
                         round(roots[1] - roots[0], 2) if roots else 0, round(float(minimum), 6)])
    return rows, [table[2] for table in tables.values()]


def write_csv(name, fields, rows):
    with (ROOT / name).open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Replace frozen tables with current Horizons results")
    args = parser.parse_args()
    rows, receipts = [], []
    for day in CASES:
        results, products = solve(day, args.refresh)
        rows.extend(results)
        receipts.extend(products)
    write_csv("contacts.csv", ["date", "spacecraft", "direction", "venus_height_km", "source_radius_in_target_radii",
                              "ingress_utc", "egress_utc", "duration_s", "minimum_limb_clearance_deg"], rows)
    write_csv("inputs.csv", ["date", "target", "epochs", "sha256_uncompressed_result_text", "query_url"], receipts)


if __name__ == "__main__":
    main()
