"""Count photons in two cleaned Swift XRT event lists for GRB 260127A.

Requires NumPy, Astropy. Run from this directory: python check_visits.py --download
"""
import argparse
import csv
import hashlib
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent
POS = (231.87925, 6.76258)  # J2000, refined XRT RA and Dec (degrees)
BASE = "https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/2026_01"
EVENTS = (
    ("01444088000", "12799867d0aeadbb074f7feefe2e86e03c8089a5a31b34e93f26756a50597409"),
    ("01444088001", "31554bde891a62c9c9fcd9f7ad950de294794b4616e950716c98fe4afbba5f76"),
)


def count_visits(path):
    with fits.open(path) as hdus:
        hdu = hdus["EVENTS"]
        hdr = hdu.header
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [hdr["TCRPX2"], hdr["TCRPX3"]]
        wcs.wcs.crval = [hdr["TCRVL2"], hdr["TCRVL3"]]
        wcs.wcs.cdelt = [hdr["TCDLT2"], hdr["TCDLT3"]]
        wcs.wcs.ctype = [hdr["TCTYP2"], hdr["TCTYP3"]]
        (x, y), = wcs.all_world2pix([POS], 1)
        events = hdu.data
        radius = np.hypot(events["X"] - x, events["Y"] - y)
        valid = (events["GRADE"] <= 12) & (events["PI"] >= 30) & (events["PI"] <= 1000)
        regions = ((10, 30), (60, 100))
        for start, stop in hdus["GTI"].data:
            active = valid & (events["TIME"] >= start) & (events["TIME"] < stop)
            counts = []
            for inner, outer in regions:
                mask = active & (radius >= inner) & (radius < outer)
                counts.extend(np.histogram(events["PI"][mask], (30, 150, 1001))[0].tolist())
            yield [hdr["OBS_ID"], round(start - hdr["TRIGTIME"], 3),
                   round(stop - hdr["TRIGTIME"], 3), *counts]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--cache", type=Path, default=HERE / "raw" / "events")
    parser.add_argument("--output", type=Path, default=HERE / "event_visits.csv")
    args = parser.parse_args()
    rows = []
    for obsid, digest in EVENTS:
        name = f"sw{obsid}xpcw3po_cl.evt.gz"
        path = args.cache / name
        if not path.exists() and args.download:
            path.parent.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(f"{BASE}/{obsid}/xrt/event/{name}", timeout=60) as response:
                raw = response.read()
            if hashlib.sha256(raw).hexdigest() != digest:
                raise ValueError(f"Archive changed: {obsid}")
            path.write_bytes(raw)
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Digest mismatch: {obsid}")
        rows.extend(count_visits(path))
    with args.output.open("w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["obs_id", "gti_start_met_minus_trigtime_s",
                         "gti_stop_met_minus_trigtime_s","source_0p3_1p5",
                         "source_1p5_10","background_0p3_1p5","background_1p5_10"])
        writer.writerows(rows)


if __name__ == "__main__":
    main()
