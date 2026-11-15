"""Check a third spacecraft at the arrival predicted by the WFC position."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np

ROOT = Path(__file__).resolve().parent
URL = "https://spdf.gsfc.nasa.gov/pub/data/ulysses/gamma/grb_1sec/1998/98351_a1.gz"
DIGEST = "b205e26e2b7533528f5ec7df286a17d2d90379e5b8dec1c1f53d77eaf13b2b13"
RECORD = np.dtype([("header", "u1"), ("subheader", "u1"), ("year_day", ">u2"),
                   ("time", ">u4"), ("count_bytes", "u1", 3), ("exposure", "u1")])


def examine(payload):
    if hashlib.sha256(payload).hexdigest() != DIGEST:
        raise ValueError("Archive file differs from the verified edition")
    rows = np.frombuffer(gzip.decompress(payload), dtype=RECORD)
    if not np.all((rows["year_day"] == 8351) & (rows["header"] == 4) & (rows["subheader"] == 0)):
        raise ValueError("Unexpected record layout or observing date")
    sec = rows["time"] / 10000
    if np.any(np.diff(sec) <= 0):
        raise ValueError("Non-monotonic times")
    counts = rows["count_bytes"].astype("i4") @ np.array([65536, 256, 1]) / 2
    live = rows["exposure"] / 8
    total = np.column_stack((counts, live)).cumsum(axis=0)
    total = np.vstack(([0.0, 0.0], total))

    def extract(a, b):
        i, j = np.searchsorted(sec, [a, b])
        return total[j] - total[i]

    def contrast(start):
        ns, es = extract(start, start + 111)
        nb, eb = extract(start - 210, start - 30) + extract(start + 141, start + 321)
        if es < 0.97 * 111 or eb < 0.97 * 360:
            return None
        return float(ns / es - nb / eb), float(es), float(eb), float(ns), float(nb)

    # Geometric ICRF vector about the geocenter. Horizons uses TDB.
    # Interpolate to reception on the spacecraft, iterating the plane-wave delay.
    block = (ROOT / "inputs/ulysses_horizons.txt").read_text().split("$$SOE")[1].split("$$EOE")[0]
    data = np.array([[float(v) for v in [row.split(",")[0], *row.split(",")[2:8]]]
                     for row in block.strip().splitlines()])
    ra, dec = np.deg2rad([294.6259, 80.5517])
    direction = np.array([np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)])
    utc = 0.079960 * 86400
    # TT–UTC was 63.184 s on this date. The sub-millisecond TDB–TT difference
    # is negligible in evaluating the slowly varying spacecraft vector.
    row = data[1]
    delay = 0.0
    for _ in range(3):
        dt = (2451164.5-row[0])*86400 + utc + 63.184 + delay
        delay = -float((row[1:4] + row[4:7]*dt) @ direction) / 299792.458
    start = utc + delay
    actual = contrast(start)
    if actual is None:
        raise ValueError("Insufficient target-window coverage")
    controls = []
    for t in range(270, 86000, 111):
        if abs(t-start) <= 650:
            continue
        result = contrast(t)
        if result is not None:
            controls.append(result[0])
    exceed = int(np.sum(np.asarray(controls) >= actual[0]))
    return {"archive": URL, "sha256": DIGEST, "delay_seconds": round(delay, 3),
            "predicted_utc_seconds_of_day": round(start, 3), "source_seconds": actual[1],
            "background_seconds": actual[2], "source_counts": actual[3], "background_counts": actual[4],
            "excess_counts_per_second": round(actual[0], 3), "controls": len(controls),
            "controls_at_least_as_large": exceed, "rank_fraction_plus_one": round((exceed+1)/(len(controls)+1), 6)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--daily-file", type=Path, help="Use the exact archived gzip file offline")
    args = parser.parse_args()
    if args.daily_file:
        raw = args.daily_file.read_bytes()
    else:
        with urlopen(URL, timeout=90) as response:
            raw = response.read()
    print(json.dumps(examine(raw), indent=2))
