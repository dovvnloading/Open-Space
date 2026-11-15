"""Replay the 2005–2035 moon–moon search. Times are geocentric UTC reception."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from datetime import datetime, timedelta
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parent
SOURCE = "https://content.cld.iop.org/journals/2632-3338/5/3/69/revision1/psjad26e9t3_mrt.txt"
ORIGIN = datetime(2005, 1, 1)
MJD_ZERO = datetime(1858, 11, 17)
AU = 149597870.7


def load_ephemeris(path):
    data = path.read_bytes() if path else urlopen(SOURCE, timeout=90).read()
    digest = json.loads((ROOT / "sources.json").read_text())["table3_sha256"]
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError("Unexpected Table 3 digest")
    rows = [line.split() for line in data.decode().splitlines() if line[:1].isdigit()]
    values = np.array([row[3:] for row in rows], dtype=float)
    times = [(datetime.fromisoformat(row[1]+" "+row[2].zfill(8))-ORIGIN).total_seconds()/86400
             for row in rows]
    t = np.array(times)
    assert np.max(abs(np.diff(t)-1/3)) < 1.1/86400  # UTC leap-second days
    return t, values


def table_minima(t, values):
    offsets = (values[:, :2]-values[:, 6:8])*1000
    interp = CubicSpline(t, offsets)
    out = []
    for i in find_peaks(-np.linalg.norm(offsets, axis=1))[0]:
        opt = minimize_scalar(lambda x: np.linalg.norm(interp(x)),
                              bounds=(t[i-1], t[i+1]), method="bounded")
        out.append((ORIGIN+timedelta(days=opt.x), float(opt.fun)))
    return out


def archive_tables():
    with gzip.open(ROOT/"archive-times.csv.gz", "rt") as stream:
        rows = list(csv.DictReader(stream))
    return {provider: [(float(r["mjd"]), r["observation"]) for r in rows
                       if r["provider"] == provider] for provider in ("CADC", "MPC")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, help="local publisher Table 3 for offline replay")
    args = parser.parse_args()
    t, values = load_ephemeris(args.table)
    minima = table_minima(t, values)
    archives = archive_tables()
    with gzip.open(ROOT/"horizons-samples.csv.gz", "rt") as stream:
        samples = list(csv.DictReader(stream))
    results = []
    for when, impact in minima:
        if impact >= 26:
            continue
        key = when.strftime("%Y-%m-%d")
        grid = [s for s in samples if s["event_day"] == key]
        # The Horizons files contain independently timed apparent directions for each moon.
        ts = np.array([float(s["mjd"]) for s in grid])
        xy = np.array([[float(s["hi_ra"]),float(s["hi_dec"]),float(s["na_ra"]),
                        float(s["na_dec"])] for s in grid])
        diff = (xy[:, 2:]-xy[:, :2])*3.6e6
        diff[:, 0] *= np.cos(np.deg2rad(xy[:, 1]))
        curve = CubicSpline(ts-ts[0], diff)
        j = np.linalg.norm(diff,axis=1).argmin()
        if j == 0 or j == len(ts)-1:
            htime, himpact = ts[j], float(np.linalg.norm(diff[j]))
            edge = True
        else:
            fit = minimize_scalar(lambda x: np.linalg.norm(curve(x)),
                                  bounds=(ts[j-1]-ts[0],ts[j+1]-ts[0]),method="bounded")
            htime, himpact, edge = ts[0]+fit.x, fit.fun, False
        mjd = (when-MJD_ZERO).total_seconds()/86400
        dist = float(grid[j]["hi_range_au"])
        result = dict(utc=when.isoformat(timespec="minutes"), nominal_mas=round(impact,3),
                      nominal_km=round(impact/206264806.247*dist*AU,1),
                      horizons_mas=round(himpact,3), horizons_km=round(himpact/206264806.247*dist*AU,1),
                      horizons_offset_hr=round(24*(htime-mjd),2), horizons_at_grid_edge=edge,
                      foreground_horizons=("Namaka" if float(grid[j]["na_range_au"])<dist else "Hiiaka"))
        for provider, records in archives.items():
            tm, name = min(records, key=lambda r:abs(r[0]-mjd))
            result[provider+"_nearest_id"] = name
            result[provider+"_offset_hr"] = round((tm-mjd)*24,2)
        results.append(result)
    with (ROOT/"events.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    future = min((b,w) for w,b in minima if w >= datetime(2026,11,15))
    print(f"{len(t)} epochs, {len(minima)} minima, {len(results)} below 26 mas;"
          f" closest future nominal separation {future[0]:.2f} mas on {future[1]:%Y-%m-%d}")


if __name__ == "__main__":
    main()
