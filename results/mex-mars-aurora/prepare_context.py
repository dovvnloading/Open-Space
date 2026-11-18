"""Reduce the archived Earth GNSS map and Mars Express radiation counters."""

import argparse
from bisect import bisect_left
import csv
from datetime import datetime
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RECEIPTS = {
    "edac": {
        "url": "https://zenodo.org/api/records/13364069/files/MEX_NDMW0D0G_MEX_NACW0D0G_2024_04_11_07_37_22.272.txt/content",
        "sha256": "9ea6f1669e132aad900ac051cc86283b64b6f7710ba9adcdea0bd3ca065a01a5",
    },
    "tec": {
        "url": "https://spdf.gsfc.nasa.gov/pub/data/gps/tec1hr_igs/2024/gps_tec1hr_igs_20240317_v01.cdf",
        "sha256": "fe1fd2fdf0ba2793a653e30f345a49a915f72ea53393c75caa1a3de0a5a523bc",
    },
}
PASSES = {
    "4432": ("2024-03-17T04:29:21", "2024-03-17T04:51:30"),
    "4433": ("2024-03-17T11:13:48", "2024-03-17T11:48:53"),
    "4434": ("2024-03-18T01:27:07", "2024-03-18T01:49:16"),
}


def verified(path, name):
    if sha256(path.read_bytes()).hexdigest() != RECEIPTS[name]["sha256"]:
        raise ValueError(f"Context source does not match the pinned digest: {path}")


def read_edac(path):
    verified(path, "edac")
    samples = {"NDMW0D0G": [], "NACW0D0G": []}
    with path.open(newline="", encoding="ascii") as stream:
        for row in csv.reader(stream, delimiter="\t"):
            if len(row) < 3 or not row[0].startswith("2024-"):
                continue
            time = datetime.fromisoformat(row[0]).isoformat(timespec="milliseconds")
            for readings, value in zip(samples.values(), row[1:]):
                if value.strip():
                    readings.append((time, int(float(value))))
    result = {}
    for name, series in samples.items():
        windows = {}
        times = [time for time, value in series]
        for volume, (start, stop) in PASSES.items():
            i, j = bisect_left(times, start), bisect_left(times, stop)
            windows[volume] = {
                "inside": series[i:j],
                "before": series[i - 1],
                "after": series[j],
            }
        changes = [(previous, current) for previous, current in zip(series, series[1:])
                   if previous[1] != current[1]]
        result[name] = {
            "passes": windows, "decreases_in_source": sum(a[1] > b[1] for a, b in changes),
            "march_17_18_changes": [(a, b) for a, b in changes if b[0][:10] in ("2024-03-17", "2024-03-18")],
        }
    return result


def read_tec(path):
    import cdflib
    import numpy as np

    verified(path, "tec")
    cdf = cdflib.CDF(path)
    latitude, longitude = 40.427204, -4.2507
    x, y = (longitude + 5) / 5, (latitude - 40) / 2.5
    weights = np.array([[(1 - y) * (1 - x), (1 - y) * x], [y * (1 - x), y * x]])
    lat = cdf.varget("lat")
    lon = cdf.varget("lon")
    indices = [int(np.flatnonzero(lat == v)[0]) for v in (40, 42.5)]
    neighbors = [int(np.flatnonzero(lon == v)[0]) for v in (-5, 0)]
    epochs = cdflib.cdfepoch.to_datetime(cdf.varget("Epoch"))
    result = {}
    for index, epoch in enumerate(epochs):
        stamp = str(epoch)[:19]
        if stamp[:13] not in {f"2024-03-17T{hour:02}" for hour in range(9, 14)}:
            continue
        result[stamp] = {}
        for key in ("tecUHR", "rmsUHR", "tecCOD", "rmsCOD", "tecCOR", "rmsCOR", "tecEHR", "rmsEHR"):
            raster = cdf.varget(key, startrec=index, endrec=index)[0]
            values = raster[np.ix_(indices, neighbors)]
            valid = np.isfinite(values) & (values >= 0) & (values <= 300)
            result[stamp][key] = {
                "neighbors": [[float(v) if keep else None for v, keep in zip(row, mask)]
                              for row, mask in zip(values, valid)],
                "bilinear_tecu": float(np.sum(weights * values)) if valid.all() else None,
            }
    return {"station_lat_lon": [latitude, longitude], "vertical_tec": result}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edac", type=Path, required=True)
    parser.add_argument("--tec", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {
        "sources": RECEIPTS,
        "edac_naive_archive_times": read_edac(args.edac),
        "earth_gnss": read_tec(args.tec),
    }
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
