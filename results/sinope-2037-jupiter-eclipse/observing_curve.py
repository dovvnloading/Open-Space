"""Rebuild the local event light curve from a compact, offline JPL state excerpt."""

import argparse
import csv
import io
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline

import search

BODIES = [10, 399, 509, 599]
STATEFILE = search.DIRECTORY / "event-states.csv"


def export(cache, path):
    with path.open("w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(
            ["id", "jd_tdb", "x_km", "y_km", "z_km", "vx_km_s", "vy_km_s", "vz_km_s"]
        )
        for body in BODIES:
            search.state_table(cache, body)
            data = (
                (cache / f"{body}.txt")
                .read_text()
                .split("$$SOE")[1]
                .split("$$EOE")[0]
                .strip()
            )
            for row in csv.reader(io.StringIO(data)):
                if 6496 <= float(row[0]) - search.J0 <= 6504:
                    writer.writerow(
                        [body, row[0].strip(), *[v.strip() for v in row[2:8]]]
                    )


def interpolants(path):
    table = np.genfromtxt(path, delimiter=",", skip_header=1)
    states = {}
    for body in BODIES:
        r = table[table[:, 0] == body]
        states[body] = CubicHermiteSpline(
            r[:, 1] - search.J0, r[:, 2:5], r[:, 5:8] * search.DAY, extrapolate=False
        )
    return states


def run(path, outdir):
    states = interpolants(path)
    mid = 6499.79801942
    with (outdir / "observing-curve.csv").open("w") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["earth_utc", "sunlight_fraction", "separation_arcmin"])
        for minutes in range(-15 * 60, 15 * 60 + 1, 10):
            t = mid + minutes / 1440
            arrival = search.earth_time(states, 509, t)
            earth = states[399](arrival)
            jupiter = search.seen_from(earth, states[599], arrival)
            sep = search.angle(states[509](t) - earth, jupiter) * 60
            w.writerow(
                [
                    search.utc(arrival),
                    round(1 - search.fraction(states, 509, 599, t), 6),
                    round(sep, 4),
                ]
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--export-from-cache",
        type=Path,
        help="Regenerate the ephemeris excerpt from the full search inputs",
    )
    parser.add_argument("--states", type=Path, default=STATEFILE)
    parser.add_argument("--outdir", type=Path, default=search.DIRECTORY)
    a = parser.parse_args()
    if a.export_from_cache:
        export(a.export_from_cache, a.states)
    run(a.states, a.outdir)
