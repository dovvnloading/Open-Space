"""Optional online discovery screen, using geometric simultaneous rather than retarded states.

Writes the cached four-year selection in .uranian-cache and a table with provisional
candidates on stdout. Use replay.py for the physical event circumstances.
"""
import csv
import io
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations, permutations
from pathlib import Path
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq, minimize_scalar

JD0 = 2469076.5  # 2048-01-01 00:00 TDB
RADII = {701: 578.9, 702: 584.7, 703: 788.9, 704: 761.4, 705: 235.8}
CACHE = Path(".uranian-cache")


def fetch(job):
    body, year = job
    path = CACHE / f"{body}_{year}.npz"
    if path.exists():
        return
    params = dict(format="json", COMMAND=f"'{body}'", CENTER="'500@7'", MAKE_EPHEM="'YES'",
                  EPHEM_TYPE="'VECTORS'", START_TIME=f"'{year}-01-01'", STOP_TIME=f"'{year+1}-01-01'",
                  STEP_SIZE="'3 h'", OUT_UNITS="'KM-S'", VEC_TABLE="'2'", CSV_FORMAT="'YES'",
                  REF_PLANE="'FRAME'", OBJ_DATA="'YES'")
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(params)
    for attempt in range(4):
        try:
            with urlopen(url, timeout=120) as reply:
                data = json.load(reply)["result"].split("$$SOE")[1].split("$$EOE")[0]
            break
        except URLError:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    arr = np.loadtxt(io.StringIO(data), delimiter=",", usecols=(0, 2, 3, 4, 5, 6, 7))
    np.savez_compressed(path, a=arr)


def interpolants():
    CACHE.mkdir(exist_ok=True)
    bodies = [*RADII, 399, 10]
    years = range(2048, 2052)
    with ThreadPoolExecutor(max_workers=5) as pool:
        list(pool.map(fetch, [(b, y) for y in years for b in bodies]))
    states = {}
    for b in bodies:
        h = np.vstack([np.load(CACHE / f"{b}_{y}.npz")["a"][:-1] for y in years])
        assert h[0, 0] == JD0
        states[b] = CubicHermiteSpline((h[:, 0] - JD0) * 86400, h[:, 1:4], h[:, 4:])
    return states


def screen_pair(states, grid, front, back, source):
    def edge(t):
        q = states[back](t)
        ray = states[source](t) - q
        dist = np.linalg.norm(ray, axis=-1)
        ray = ray / dist[..., None]
        delta = states[front](t) - q
        depth = np.sum(delta * ray, axis=-1)
        sep = np.linalg.norm(np.cross(delta, ray), axis=-1)
        width = RADII[front] + RADII[back] + (695700 * depth / dist if source == 10 else 0.0)
        return sep - width, depth

    value, depth = edge(grid)
    minima = 1 + np.flatnonzero((depth[1:-1] > 0) & (value[1:-1] < value[:-2]) &
                                (value[1:-1] <= value[2:]) & (value[1:-1] < 3500))
    events = []
    for k in minima:
        def residual(t):
            return edge(t)[0]
        mid = minimize_scalar(residual, bounds=(grid[k-1], grid[k+1]), method="bounded",
                              options={"xatol": 0.01}).x
        if residual(mid) >= 0:
            continue
        def contact(sign):
            width = 1200
            while residual(mid + sign * width) < 0:
                width *= 2
            return brentq(residual, *sorted([mid, mid + sign * width]))
        events.append((front, back, source, contact(-1), contact(1)))
    return events


def main():
    d = interpolants()
    # Last 3 hours are omitted because year-seam duplicates were dropped.
    grid = np.arange(0, 1461*86400 - 10800, 600.)
    events = [event for source in (10, 399) for front, back in permutations(RADII, 2)
              for event in screen_pair(d, grid, front, back, source)]
    pairs = [(a, b) for a, b in combinations(events, 2)
             if a[:2] == b[:2] and a[2] != b[2] and max(a[3], b[3]) < min(a[4], b[4])]
    print(f"{len(events)} intervals. {len(pairs)} simultaneous same-pair candidates.", file=sys.stderr)
    writer = csv.writer(sys.stdout)
    writer.writerow(["foreground", "background", "overlap_start_jd_tdb", "overlap_end_jd_tdb"])
    for a, b in pairs:
        writer.writerow([a[0], a[1], JD0 + max(a[3], b[3])/86400, JD0 + min(a[4], b[4])/86400])


if __name__ == "__main__":
    main()
