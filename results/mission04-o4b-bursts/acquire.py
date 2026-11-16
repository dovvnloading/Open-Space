"""Prepare fixed O4b H/L survey epochs and validated GWOSC 4 kHz R1 data.

Run: python acquire.py PATH_TO_CACHE
Only DATA coverage and civil month boundaries determine survey file choices.
"""
from __future__ import annotations

import argparse
import calendar
import hashlib
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

BASE = "https://gwosc.org"
RELEASE = "O4b_4KHZ_R1"
DURATION = 4096
FS = 4096
RUN_END = 1422118818
DETECTORS = ("H1", "L1", "V1")
CONTROL = "GW250114_082203-v2"
DQ_REQUIRED = ("DATA", "BURST_CAT1", "BURST_CAT2", "BURST_CAT3")
INJ_REQUIRED = ("NO_CBC_HW_INJ", "NO_BURST_HW_INJ", "NO_DETCHAR_HW_INJ", "NO_STOCH_HW_INJ")
META = {
    "MD5.txt": f"{BASE}/archive/md5/{RELEASE}/strain-hdf.txt",
    "allevents.json": f"{BASE}/eventapi/json/allevents/",
    **{f"{det}_DATA.txt": f"{BASE}/timeline/segments/O4b_16KHZ_R1/{det}_DATA/1396417050/25701768//"
       for det in DETECTORS},
}


def fetch(url, path, checksum=None):
    """Cache downloads atomically. Always verify archive MD5 when supplied."""
    def valid(file):
        if not file.exists():
            return False
        if checksum is None:
            return True
        with file.open("rb") as stream:
            return hashlib.file_digest(stream, "md5").hexdigest() == checksum

    if valid(path):
        return
    tmp = path.with_name(path.name + ".part")
    try:
        subprocess.run(["curl", "--fail", "--silent", "--show-error", "--location",
                        "--retry", "3", "--retry-all-errors", "--max-time", "240",
                        url, "--output", str(tmp)], check=True)
        if not valid(tmp):
            raise ValueError(f"GWOSC checksum mismatch for {path.name}")
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def gps_midnight(year, month, day):
    # GPS minus UTC is 18 s throughout O4b.
    return int(datetime(year, month, day, tzinfo=timezone.utc).timestamp()) - 315964800 + 18


def select_blocks(cache):
    segments = {det: np.loadtxt(cache / f"{det}_DATA.txt", dtype=int, ndmin=2)[:, :2]
                for det in DETECTORS}
    records = {}
    for line in (cache / "MD5.txt").read_text().splitlines():
        checksum, path = line.split()
        detector, directory, name = path.split("/")
        start = int(name.split("-")[-2])
        records[(detector, start)] = {"md5": checksum, "url": f"{BASE}/archive/data/{RELEASE}/{directory}/{name}"}

    events = json.loads((cache / "allevents.json").read_text())["events"]
    control_gps = events[CONTROL]["GPS"]
    blocks = []
    for n in range(10):
        year, month = divmod(2024 * 12 + 3 + n, 12)
        month += 1
        low = gps_midnight(year, month, 15)
        high = min(RUN_END, gps_midnight(year, month, calendar.monthrange(year, month)[1]) + 86400)
        first = (low + DURATION - 1) // DURATION * DURATION
        late = range(first, high - DURATION + 1, DURATION)
        month_start = max(1396796418, gps_midnight(year, month, 1))
        early = range((low - DURATION) // DURATION * DURATION, month_start - 1, -DURATION)
        # If no later block exists in the month (July maintenance), take the
        # last complete block ending before the 15th. This uses only DATA.
        for gps in [*late, *early]:
            if all(np.any((segments[d][:, 0] <= gps) & (segments[d][:, 1] >= gps + DURATION))
                   for d in ("H1", "L1")):
                blocks.append({"role": "search", "month": f"{year:04d}-{month:02d}", "gps": gps,
                               "selection": "first_on_after_15th" if gps >= low else "last_before_15th"})
                break
        else:
            raise ValueError(f"No complete H/L DATA block in {year}-{month}")
    blocks.append({"role": "control", "month": "2025-01", "gps": int(control_gps) // DURATION * DURATION,
                   "event": CONTROL, "event_gps": control_gps})
    manifest = {}
    for block in blocks:
        gps = block["gps"]
        block["files"] = {}
        for d in DETECTORS:
            if (d, gps) not in records:
                if d != "V1":
                    raise ValueError(f"Missing required detector {d} at {gps}")
                continue
            record = records[(d, gps)]
            name = record["url"].rsplit("/", 1)[1]
            manifest[name] = record
            block["files"][d] = name
        block["events"] = sorted([{"name": k, "gps": v["GPS"], "catalog": v["catalog.shortName"]}
                                  for k, v in events.items() if gps-64 <= v["GPS"] < gps+DURATION+64],
                                 key=lambda x: x["gps"])
    return {"release": RELEASE, "fs": FS, "duration": DURATION, "blocks": blocks, "manifest": manifest,
            "metadata": {k: {"url": v, "sha256": hashlib.sha256((cache/k).read_bytes()).hexdigest()}
                         for k,v in META.items()}}


def require_flags(group, mask_name, short_name, required, gps):
    flags = [s.decode("ascii") for s in group[short_name][:]]
    mask = group[mask_name]
    if mask.shape != (DURATION,) or mask.attrs["Xstart"] != gps or mask.attrs["Xspacing"] != 1:
        raise ValueError(f"Invalid 1 Hz mask metadata: {mask.name}")
    bit_mask = sum(1 << flags.index(s) for s in required)
    return mask[:] & bit_mask == bit_mask


def load_detector(path):
    """Return integer GPS, raw strain and per-second DQ/injection/finite-data validity.

    Keep low-amplitude persistent CW injections in the transient search. The scan's
    fixed line treatment and candidate follow-up must account for these narrow lines.
    """
    path = Path(path)
    with h5py.File(path, "r") as f:
        gps = int(f["meta/GPSstart"][()])
        d = f["meta/Detector"][()].decode("ascii")
        expected = f"{d[0]}-{d}_GWOSC_{RELEASE}-{gps}-{DURATION}.hdf5"
        if path.name != expected or f["meta/Duration"][()] != DURATION:
            raise ValueError(f"File identity mismatch: {path.name}")
        h = f["strain/Strain"]
        if h.shape != (DURATION * FS,) or h.attrs["Xspacing"] != 1/FS or h.attrs["Xstart"] != gps:
            raise ValueError("Invalid GWOSC time grid")
        good = require_flags(f["quality/simple"], "DQmask", "DQShortnames", DQ_REQUIRED, gps)
        good &= require_flags(f["quality/injections"], "Injmask", "InjShortnames", INJ_REQUIRED, gps)
        strain = h[:]
    good &= np.isfinite(strain.reshape(DURATION, FS)).all(axis=1)
    return gps, strain, good


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    cache = args.cache
    cache.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(lambda kv: fetch(kv[1], cache/kv[0]), META.items()))
        selection = select_blocks(cache)
        (cache/"search_selection.json").write_text(json.dumps(selection, indent=2)+"\n")
        for block in selection["blocks"]:
            print(block["role"], block["month"], block["gps"], ",".join(block["files"]), flush=True)
        if args.metadata_only:
            return
        def download(item):
            name, rec = item
            fetch(rec["url"], cache/name, rec["md5"])
            print(f"Verified {name} {rec['md5']}", flush=True)
        list(pool.map(download, selection["manifest"].items()))


if __name__ == "__main__":
    main()
