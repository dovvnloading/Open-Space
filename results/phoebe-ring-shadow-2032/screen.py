"""Find Phoebe's near crossings of Saturn's antisolar axis using JPL day-grid state vectors."""
import argparse
from datetime import datetime, timedelta
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar


def states(body, cache):
    path = cache / f"saturn-centered-{body}.txt"
    if not path.exists():
        params = {"COMMAND": str(body), "CENTER": "500@699", "START_TIME": "2026-01-01",
                  "STOP_TIME": "2041-01-01", "STEP_SIZE": "1 d", "OBJ_DATA": "NO",
                  "MAKE_EPHEM": "YES", "EPHEM_TYPE": "VECTORS", "VEC_TABLE": "2",
                  "VEC_CORR": "NONE", "CSV_FORMAT": "YES", "OUT_UNITS": "KM-S", "REF_PLANE": "FRAME"}
        query = urlencode(dict(format="text", **{k: f"'{v}'" for k, v in params.items()}))
        with urlopen("https://ssd.jpl.nasa.gov/api/horizons.api?" + query, timeout=90) as response:
            raw = response.read()
        if b"$$EOE" not in raw:
            raise ValueError(raw[:300])
        path.write_bytes(raw)
    records = path.read_text().split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    rows = np.array([[float(row[k]) for k in (0, 2, 3, 4, 5, 6, 7)]
                     for row in (line.split(",") for line in records)])
    return CubicHermiteSpline(rows[:, 0], rows[:, 1:4], rows[:, 4:] * 86400)


def search(cache):
    cache.mkdir(parents=True, exist_ok=True)
    moon, sun = states(609, cache), states(10, cache)

    def impact(t):
        v, s = moon(t), sun(t)
        s = s / np.linalg.norm(s, axis=-1, keepdims=True)
        axial = np.sum(v * s, axis=-1)
        return np.linalg.norm(v - axial[..., None] * s, axis=-1), axial

    t = np.arange(moon.x[0], moon.x[-1], .125)
    dist, axial = impact(t)
    ix = np.flatnonzero((dist[1:-1] < dist[:-2]) & (dist[1:-1] < dist[2:])
                        & (dist[1:-1] < 150000) & (axial[1:-1] < 0)) + 1
    for i in ix:
        best = minimize_scalar(lambda offset: impact(t[i] + offset)[0],
                               bounds=(-.15, .15), method="bounded")
        epoch = t[i] + best.x
        date = datetime(2000, 1, 1, 12) + timedelta(days=epoch - 2451545)
        b, z = impact(epoch)
        print(date.isoformat(timespec="minutes"), f"{b:.0f} km off axis; {-z:.0f} km behind Saturn")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache", type=Path, help="local directory for two downloaded Horizons responses")
    search(p.parse_args().cache)
