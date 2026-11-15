"""Audit the archived results offline, or repeat the window scan on HEASARC FITS products."""

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
BASE = "https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/"


def load():
    with (ROOT / "observations.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    with (ROOT / "windows.csv").open() as stream:
        windows = list(csv.DictReader(stream))
    return rows, windows


def summary(rows, windows):
    histogram = Counter(int(r["max_hits"]) for r in rows)
    groups = {}
    for r in rows:
        groups.setdefault(r["field"], []).append(float(r["gti_seconds"]))
    if len({r["obsid"] for r in rows}) != len(rows):
        raise ValueError("Duplicate observation in frozen population")
    by_obs = {r["obsid"]: r for r in rows}
    for w in windows:
        r = by_obs[w["obsid"]]
        if not 5 <= int(w["hits"]) <= int(r["max_hits"]) or float(w["span_us"]) > 10.04:
            raise ValueError(f"Invalid interval in {r['obsid']}")
    if {r["obsid"] for r in rows if int(r["max_hits"]) >= 5} != {w["obsid"] for w in windows}:
        raise ValueError("Missing or unexpected threshold crossing")
    return {
        "observations": len(rows),
        "seconds": round(math.fsum(float(r["gti_seconds"]) for r in rows), 3),
        "in_band_events": sum(int(r["events"]) for r in rows),
        "observation_maxima": dict(sorted(histogram.items())),
        "threshold_windows": len(windows),
        "all_slow_windows": sum(int(w["fast_hits"]) == 0 for w in windows),
        "fields": {k: [len(v), round(math.fsum(v), 3)] for k, v in sorted(groups.items())},
    }


def scan(row, cache):
    import numpy as np
    from astropy.io import fits

    obs = row["obsid"]
    filename = f"ni{obs}_0mpu7_cl.evt.gz"
    path = cache / filename if cache is not None else None
    if path is not None and path.exists():
        data = path.read_bytes()
    else:
        url = f"{BASE}{row['archive_month']}/{obs}/xti/event_cl/{filename}"
        with urlopen(url, timeout=120) as response:
            data = response.read()
        if path is not None:
            path.write_bytes(data)
    if hashlib.sha256(data).hexdigest() != row["sha256"]:
        raise ValueError(f"Archive content changed for {obs}")
    with fits.open(io.BytesIO(gzip.decompress(data))) as hdus:
        events = hdus["EVENTS"].data
        times = events["TIME"]
        mask = (events["PI"] >= 50) & (events["PI"] <= 800)
        ix = np.flatnonzero(mask)
        ix = ix[np.argsort(times[ix], kind="stable")]
        t = times[ix]
        # Products are the pipeline-cleaned event tables. No extra background cut.
        # The published event times are binary64 MET. Adding 10 microseconds
        # rounds to the nearest stored endpoint (at most 0.030 microsecond here).
        end = np.searchsorted(t, t + 1e-5, side="right")
        counts = end - np.arange(len(t))
        gtis = hdus["GTI"].data
        exposure = math.fsum(gtis["STOP"] - gtis["START"])
        measured = (len(t), int(counts.max(initial=0)))
        expected = (int(row["events"]), int(row["max_hits"]))
        if measured != expected or abs(exposure - float(row["gti_seconds"])) > 1e-4:
            raise ValueError(f"Recomputed exposure or counts changed for {obs}")
        windows = []
        last_end = 0
        for start in np.flatnonzero(counts >= 5):
            if start < last_end:
                continue
            last_end = end[start]
            ids = ix[start:last_end]
            if not np.any((gtis["START"] <= t[start]) & (gtis["STOP"] >= t[last_end - 1])):
                raise ValueError(f"Window crosses a GTI boundary for {obs}")
            flags = events["EVENT_FLAGS"][ids]
            fast = sum(bool(flag[4]) for flag in flags)  # FITS 8X, fast is bit 3
            if fast != int(np.sum(events["PHA_FAST"][ids] >= 0)):
                raise ValueError(f"Fast chain / pulse height mismatch for {obs}")
            windows.append({
                "obsid": obs,
                "start_met": f"{t[start]:.8f}",
                "span_us": f"{1e6 * (t[last_end-1] - t[start]):.3f}",
                "hits": str(len(ids)),
                "distinct_detectors": str(len(set(events["DET_ID"][ids]))),
                "fast_hits": str(fast),
                "detector_ids": " ".join(map(str, events["DET_ID"][ids])),
                "pi_channels": " ".join(map(str, events["PI"][ids])),
            })
    return windows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--obsid", help="Re-read one observation")
    selection.add_argument("--threshold-files", action="store_true", help="Re-read all files containing a five-hit window")
    selection.add_argument("--all", action="store_true", help="Re-read every event file")
    parser.add_argument("--cache", type=Path, help="Optional directory for original compressed FITS files")
    args = parser.parse_args()
    rows, windows = load()
    print(json.dumps(summary(rows, windows), indent=2))
    if args.cache is not None:
        args.cache.mkdir(parents=True, exist_ok=True)
    if args.obsid and args.obsid not in {r["obsid"] for r in rows}:
        parser.error("ObsID absent from snapshot")
    for row in rows:
        obs = row["obsid"]
        if args.all or args.obsid == obs or (args.threshold_files and int(row["max_hits"]) >= 5):
            measured = scan(row, args.cache)
            if measured != [w for w in windows if w["obsid"] == obs]:
                raise ValueError(f"Event-level reproduction failed for {obs}")
            print(f"{obs} verified")


if __name__ == "__main__":
    main()
