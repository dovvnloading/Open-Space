"""Recheck the exact earlier QLP products that yielded usable confirmation-target coverage."""
import argparse
import gzip
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from acquire import fetch
from pixels import read_stamp
from scan import Binning, Curve, detect


def audit(job):
    tic, visit = job
    raw = fetch(visit["url"])
    if hashlib.sha256(raw).hexdigest() != visit["sha256"]:
        raise ValueError(f"Changed archive product for TIC {tic}, sector {visit['sector']}")
    h, _, table = read_stamp(raw)
    if int(h["TICID"]) != int(tic) or int(h["SECTOR"]) != int(visit["sector"]):
        raise ValueError("Archive identification mismatch")
    curve = Curve(table["TIME"].astype(float), table[visit["flux_column"]].astype(float),
                  table["QUALITY"] == 0, {}, {})
    bins = Binning(curve)
    events = detect(bins.bin_flux(bins.flux), bins.centers, require_shoulders=True)
    if len(events) != len(visit["events"]) or int(bins.valid.sum()) != visit["valid_bins"]:
        raise ValueError(f"Coverage or event count changed for TIC {tic}")
    for found, saved in zip(events, visit["events"]):
        if not np.allclose([found["peak_btjd"], found["depth"]], [saved["peak_btjd"], saved["depth"]],
                           rtol=0, atol=1e-7):
            raise ValueError(f"Event measurement changed in TIC {tic}")
    return len(events)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("receipts", type=Path)
    args = ap.parse_args()
    with gzip.open(args.receipts, "rt") as stream:
        stars = json.load(stream)
    jobs = [(s["tic"], v) for s in stars for v in s["visits"] if v["status"] == "ok"]
    with ThreadPoolExecutor(5) as pool:
        counts = list(pool.map(audit, jobs))
    print(f"{len(jobs)} previous visits verified, with {sum(counts)} screening events")
