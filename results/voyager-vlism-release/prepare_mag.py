"""Extract frozen PWS rows and matched definitive V1 MAG windows."""

import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path

import cdflib
import numpy as np

from audit_catalog import ROOT, archive_members, format_utc, product_rows, utc


START = utc("2023-01-01T00:00:00Z")
STOP = utc("2025-05-01T00:00:00Z")
UNIX_EPOCH = datetime(1970, 1, 1, tzinfo=timezone.utc)


def positive(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) and number > 0 else None


def read_mag(directory):
    sources = json.loads((ROOT / "mag-sources.json").read_text(encoding="utf-8"))
    times, fields = [], []
    for source in sources:
        path = directory / source["file"]
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        if path.stat().st_size != source["bytes"] or digest.hexdigest() != source["raw_sha256"]:
            raise ValueError(f"MAG source mismatch: {path}")
        cdf = cdflib.CDF(path)
        epoch = cdf.varget("Epoch")
        f1 = np.asarray(cdf.varget("F1"), dtype=float)
        attrs = cdf.varattsget("F1")
        valid = np.isfinite(f1) & (f1 > 0) & (f1 < 10) & (f1 != attrs["FILLVAL"])
        # Epoch labels the start of the 48-second integration.
        midpoint = cdflib.cdfepoch.to_datetime(epoch).astype("datetime64[ms]") + np.timedelta64(24, "s")
        valid &= ~np.isnat(midpoint)
        times.append(midpoint[valid].astype(np.int64))
        fields.append(f1[valid])
    time = np.concatenate(times)
    field = np.concatenate(fields)
    order = np.argsort(time, kind="stable")
    time, field = time[order], field[order]
    if np.any(np.diff(time) <= 0):
        raise ValueError("Duplicate or unordered MAG interval midpoints")
    return sources, time, field


def window(time, field, center, hours, minimum):
    half = hours * 60 * 60 * 1000
    left = np.searchsorted(time, center - half, side="left")
    right = np.searchsorted(time, center + half, side="left")
    count = int(right - left)
    return {"count": count, "field_nt": float(np.median(field[left:right])) if count >= minimum else None}


def run(directory):
    mag_files, time, field = read_mag(directory)
    records = product_rows(archive_members(), 1)
    output = []
    for row in records:
        instant = utc(row["SCET"])
        if not START <= instant < STOP:
            continue
        center = (instant - UNIX_EPOCH) // timedelta(milliseconds=1)
        density = positive(row["N_e"])
        low, high = positive(row["N_e_min"]), positive(row["N_e_max"])
        if density is None or low is None or high is None or not low <= density <= high:
            low = high = None
        output.append({
            "scet": format_utc(instant), "source": row["Source"].lower(),
            "density": density, "low": low, "high": high,
            "embedded": positive(row["Bmag"]),
            "w12": {str(day): window(time, field, center + day * 86400000, 12, 450)
                    for day in (0, -14, -7, 7, 14)},
            "w6": {str(day): window(time, field, center + day * 86400000, 6, 225)
                   for day in (0, -14, -7, 7, 14)},
            "embedded_w12": {str(day): window(time, field, center + int(day * 86400000), 12, 450)
                             for day in (-0.5, 0, 0.5)},
        })
    return {"mag_files": mag_files, "rows": output}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(args.data_root), indent=2, allow_nan=False) + "\n",
                           encoding="utf-8", newline="\n")
