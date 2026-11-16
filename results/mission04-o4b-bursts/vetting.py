"""Fixed O4b public-event exclusions and provenance. Times are GPS seconds."""
import argparse
import csv
from functools import lru_cache
import hashlib
import io
import json
import math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np


START, STOP = 1396796418, 1422118818
MARGIN = 32
HERE = Path(__file__).resolve().parent
FILENAME = "public_events.csv"
RECEIPT_FILENAME = "public_event_receipts.json"
GWTC_URL = ("https://zenodo.org/records/20348004/files/"
            "IGWN-GWTC5p0-29ebe06b7_25-SearchSummaryTable.hdf5?download=1")
GWTC_MD5 = "d4c0e5e32bb00139eaaaa15fb810d95b"
GWOSC_URL = "https://gwosc.org/eventapi/json/allevents/"
GRACE_URL = "https://gracedb.ligo.org/api/v2/superevents/?" + urlencode(
    {"query": f"gpstime: {START} .. {STOP}", "count": 1000}
)
COLUMNS = ["gps", "source", "id", "detail"]


def sha256(data):
    return hashlib.sha256(data).hexdigest()


@lru_cache(maxsize=1)
def events():
    """Use the frozen event list even when remote catalogs later change."""
    raw = (HERE / FILENAME).read_bytes()
    receipt = json.loads((HERE / RECEIPT_FILENAME).read_text())
    if sha256(raw) != receipt["csv_sha256"]:
        raise ValueError("Pinned O4b public event table checksum mismatch")
    return tuple({**row, "gps": float(row["gps"])}
                 for row in csv.DictReader(io.StringIO(raw.decode("utf-8"))))


def matches(start, end):
    """Rows with trigger centers in the supplied half-open GPS interval."""
    return [row for row in events() if start <= row["gps"] < end]


def mask_seconds(gps_start, valid_1hz):
    """Return a new 1 Hz mask rejecting seconds intersecting [t-32,t+32)."""
    valid = np.array(valid_1hz, dtype=bool, copy=True)
    if valid.ndim != 1 or gps_start != int(gps_start):
        raise ValueError("Expected a 1D one-Hz mask aligned to integer GPS")
    stop = gps_start + len(valid)
    for event in matches(gps_start - MARGIN, stop + MARGIN):
        offset = event["gps"] - gps_start
        lo = max(0, math.floor(offset - MARGIN))
        hi = min(len(valid), math.ceil(offset + MARGIN))
        valid[lo:hi] = False
    return valid


def rebuild(cache, output):
    """Read cached bytes or download once. Preserve source hashes and all versions."""
    import h5py

    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    inputs, rows = [], []

    def acquire(url, name):
        path = cache / name
        if not path.exists():
            with urlopen(url, timeout=120) as response:
                raw = response.read()
            path.write_bytes(raw)
        else:
            raw = path.read_bytes()
        inputs.append({"url": url, "bytes": len(raw), "sha256": sha256(raw)})
        return raw

    def add(t, source, identifier, detail):
        if START <= t < STOP:
            rows.append((t, source, identifier, detail))

    raw = acquire(GWTC_URL, "gwtc.hdf5")
    if hashlib.md5(raw).hexdigest() != GWTC_MD5:
        raise ValueError("Inconsistent with the LVK-published HDF5 checksum")
    with h5py.File(io.BytesIO(raw), "r") as handle:
        table = handle["search_summary"][:]
    if len(table) != 1886:
        raise ValueError("Unexpected GWTC-5.0 public candidate count")
    for row in table:
        add(float(row["gps_time"]), "LVK", row["gw_name"].decode(),
            row["dq_status"].decode())

    catalog = json.loads(acquire(GWOSC_URL, "allevents.json"))["events"]
    for identifier, row in catalog.items():
        add(float(row["GPS"]), "GWOSC", identifier, row["catalog.shortName"])

    url, page, total = GRACE_URL, 0, None
    ids = set()
    while url:
        result = json.loads(acquire(url, f"page-{page}.json"))
        if total is None:
            total = result["numRows"]
        if result["numRows"] != total:
            raise ValueError("GraceDB count changed during pagination")
        for row in result["superevents"]:
            identifier = row["superevent_id"]
            if identifier in ids or row["category"].lower() != "production":
                raise ValueError("Unexpected GraceDB result")
            ids.add(identifier)
            flags = set(row["labels"])
            detail = ",".join(sorted(flags & {
                "DQOK", "DQV", "EM_Selected", "ADVREQ", "ADVOK", "ADVNO",
                "RETRACTED", "INJ", "LOW_SIGNIF_LOCKED", "SIGNIF_LOCKED",
            }))
            add(float(row["t_0"]), "GraceDB", identifier, detail)
        page += 1
        url = result["links"]["next"]
    if len(ids) != total or total < 1:
        raise ValueError("Incomplete GraceDB pagination")

    csv_path = output / FILENAME
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(sorted(rows))
    receipt = {
        "bounds_gps": [START, STOP], "margin_seconds": MARGIN,
        "full_records": {"LVK": len(table), "GWOSC": len(catalog), "GraceDB": total},
        "rows_in_bounds": {s: sum(r[1] == s for r in rows) for s in ("LVK", "GWOSC", "GraceDB")},
        "csv_sha256": sha256(csv_path.read_bytes()), "inputs": inputs,
    }
    (output / RECEIPT_FILENAME).write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rebuild(args.cache, args.output)
