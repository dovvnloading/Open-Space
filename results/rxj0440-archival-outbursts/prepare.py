"""Rebuild frozen CSVs using source bytes from the manifest (NumPy, Astropy for FITS)."""
import csv
import gzip
import hashlib
import io
import json
import sys
from pathlib import Path
from urllib.request import urlopen

import numpy as np

ROOT = Path(__file__).resolve().parent


def csv_bytes(header, rows):
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return output.getvalue().encode("ascii")


def pack(name, payload):
    output = ROOT / "data" / name
    if output.suffix == ".gz":
        with output.open("wb") as f, gzip.GzipFile(filename="", fileobj=f, mode="wb", mtime=0) as g:
            g.write(payload)
    else:
        output.write_bytes(payload)


def prepare(paths, table):
    for k, columns, header in [
        ("bat", [0, 1, 2, 7], ["mjd", "rate", "error", "data_flag"]),
        ("maxi", list(range(9)), ["mjd", "flux", "error", "flux_2_4", "error_2_4",
                                 "flux_4_10", "error_4_10", "flux_10_20", "error_10_20"]),
    ]:
        a = np.loadtxt(paths[k])[:, columns]
        pack(k + "-daily.csv.gz", csv_bytes(header, a))
        subset = a[(a[:, 0] >= 55207) & (a[:, 0] <= 55267)]
        pack(k + "-february.csv", csv_bytes(header, subset))
    header = ["range_start_mjd", "range_stop_mjd", "barytime_mjd", "detected",
              "frequency_hz", "frequency_error_hz", "amplitude", "amplitude_error"]
    rows = ([*r["TRANGE"], r["BARYTIME"], int(r["DETECTED"]), r["FREQUENCY"],
             r["FREQUENCY_ERR"], float(r["AMPLITUDE"]), float(r["AMPSIG"])] for r in table)
    pack("gbm-detections.csv", csv_bytes(header, rows))


if __name__ == "__main__":
    from astropy.io import fits

    cache = Path(sys.argv[1])
    cache.mkdir(exist_ok=True, parents=True)
    paths = {}
    for name, source in json.loads((ROOT / "sources.json").read_text()).items():
        path = cache / source["url"].split("/")[-1]
        if not path.exists():
            with urlopen(source["url"], timeout=90) as inp:
                path.write_bytes(inp.read())
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != source["sha256"]:
            raise ValueError(f"Archive snapshot mismatch for {name}: {digest}")
        paths[name] = path
    with fits.open(paths["gbm"]) as hdus:
        prepare(paths, hdus["PSR_HIST"].data)
