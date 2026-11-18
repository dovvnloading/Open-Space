"""Check native Juno products and extract the fields needed for the frozen test."""

import argparse
import csv
from datetime import datetime
import hashlib
import io
import json
from pathlib import Path
import re

import numpy as np
from scipy.optimize import minimize_scalar
import spiceypy as spice

ROOT = Path(__file__).resolve().parent
REFERENCE = {81: "2026-02-25T09:30:00", 82: "2026-03-29T22:30:00", 83: "2026-05-01T15:00:00"}
PIXELS = (1, 3, 2, 4, 5)
WIDTH = 9622


def verified(path, digest, size=None):
    data = path.read_bytes()
    if (size is not None and len(data) != size) or hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"Wrong archive bytes at {path}")
    return data


def load_sources(root):
    receipt = json.loads((ROOT / "sources.json").read_text(encoding="utf-8"))
    for item in receipt:
        verified(root / item["path"], item["sha256"], item["bytes"])
        if item.get("spice"):
            spice.furnsh(str(root / item["path"]))
    return receipt


def column_slices(format_text):
    fields = {}
    for block in re.findall(r"^\s*OBJECT\s*=\s*COLUMN\s*$(.*?)^\s*END_OBJECT",
                            format_text, re.M | re.S):
        name = re.search(r'NAME\s*=\s*"([^"]+)"', block).group(1)
        start = int(re.search(r"START_BYTE\s*=\s*(\d+)", block).group(1)) - 1
        width = int(re.search(r"BYTES\s*=\s*(\d+)", block).group(1))
        fields[name] = slice(start, start + width)
    return fields


def read_radiation(data, fields):
    headers = [fields[name] for name in
               ("ETLookDirection", "Duration", "SectorDeadTime", "LargePixel", "Quality", "Sector")]
    channels = [[fields[f"T{pixel}IXF{channel:02} COUNTS"] for channel in range(2, 24)]
                for pixel in PIXELS]
    for line in data.splitlines(keepends=True):
        if len(line) != WIDTH:
            raise ValueError(f"Unexpected JEDI fixed record width of {len(line)}")
        time, duration, dead, large, quality, sector = [line[span] for span in headers]
        mode = large.strip().lower()
        if mode not in (b"true", b"false"):
            raise ValueError(f"Unexpected detector mode {mode!r}")
        rates = []
        for pixel in channels:
            counts = [float(line[span]) for span in pixel]
            rates.append(sum(counts) if all(np.isfinite(x) and x >= 0 for x in counts) else np.nan)
        yield (float(time), float(duration), float(dead), mode == b"true",
               int(quality), int(sector), *rates)


def decode_jedi(root, day):
    base = f"JED_180_HIERSISP_CDR_2026{day}_V04"
    directory = root / f"{base}.ranges"
    receipt = json.loads((ROOT / f"ranges-{day}.json").read_text(encoding="utf-8"))
    fields = column_slices((root / "jedi-range-metadata" / f"FMT_{day}.txt").read_text())
    records = []
    next_start = receipt["first_row"]
    for item in receipt["ranges"]:
        start, stop = item["first_row"], item["stop_row"]
        if start != next_start or stop <= start:
            raise ValueError("JEDI ranges must cover consecutive rows")
        contents = verified(directory / item["part"], item["sha256"], (stop - start) * WIDTH)
        rows = list(read_radiation(contents, fields))
        if len(rows) != stop - start:
            raise ValueError("JEDI range row count mismatch")
        records.extend(rows)
        next_start = stop
    if next_start != receipt["stop_row"]:
        raise ValueError("JEDI final range stops early")
    return np.array(records)


def closest(target, reference, hours):
    center = spice.str2et(reference)
    def radius(offset):
        return np.linalg.norm(spice.spkpos("-61", center + offset, "J2000", "NONE", target)[0])
    offsets = np.arange(-3600 * hours, 3600 * hours + 1, 60)
    best = int(np.argmin([radius(offset) for offset in offsets]))
    if best in (0, len(offsets) - 1):
        raise ValueError("Minimum must be bracketed")
    fit = minimize_scalar(radius, bounds=(offsets[best - 1], offsets[best + 1]),
                          method="bounded", options={"xatol": 0.001})
    return {"et": center + fit.x, "distance_km": fit.fun,
            "utc": spice.et2utc(center + fit.x, "ISOC", 3)}


def decode_sru(root, orbit, dates):
    table = root / f"SRU_CRT_{orbit}_{dates}_V01.CSV"
    label = table.with_suffix(".LBL")
    declared = int(re.search(r"^\s*ROWS\s*=\s*(\d+)", label.read_text(), re.M).group(1))
    rows = []
    for row in csv.reader(io.StringIO(table.read_text())):
        if len(row) != 6:
            raise ValueError(f"Incorrect SRU field count in {table}")
        epoch, clock, sensor, cluster, dn, derived = (value.strip() for value in row)
        if not epoch[0].isdigit():
            raise ValueError("Unexpected SRU record")
        cluster, dn, derived = (float(value) for value in (cluster, dn, derived))
        cluster_ok, dn_ok = cluster > 1e-30, dn > 1e-30
        method = "cluster" if cluster_ok and not dn_ok else "dn" if dn_ok and not cluster_ok else None
        rows.append({"utc": epoch, "et": spice.str2et(epoch), "sensor": int(sensor),
                     "method": method, "derived": derived if derived > 1e-30 else None})
    return {"declared_rows": declared, "images": rows}


def decode_asc(root, orbit, date):
    table = root / f"ASC_RAD_L3_PJ{orbit:03}_{date}_CHUC_V01.TAB"
    rows = []
    for line in table.read_bytes().splitlines():
        counts, bbo, seconds = (float(line[start:stop]) for start, stop in ((16,21),(23,24),(43,57)))
        rows.append((seconds, counts, bbo))
    return np.asarray(rows)


def run(root, bundle, geometry_path):
    load_sources(root)
    geometry = {}
    for orbit, dates in ((81,"2026055_2026056"), (82,"2026088_2026089"), (83,"2026120_2026121")):
        geometry[str(orbit)] = {
            "jupiter": closest("599", REFERENCE[orbit], 2),
            "thebe": closest("514", REFERENCE[orbit], 10),
            "sru": decode_sru(root, orbit, dates),
        }
    for item in geometry.values():
        for row in item["sru"]["images"]:
            time = datetime.fromisoformat(row["utc"])
            row["utc_seconds_since_j2000_noon"] = (
                time - datetime(2000, 1, 1, 12)
            ).total_seconds()
    np.savez_compressed(
        bundle, jedi_82=decode_jedi(root, "088"), jedi_83=decode_jedi(root, "121"),
        asc_81=decode_asc(root, 81, "2026056T093254"),
        asc_82=decode_asc(root, 82, "2026088T224714"),
    )
    geometry_path.write_text(json.dumps(geometry, allow_nan=False, indent=2) + "\n",
                             encoding="utf-8", newline="\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--geometry", type=Path, required=True)
    args = parser.parse_args()
    run(args.data_root, args.output, args.geometry)
