"""Freeze the public JPL Horizons vectors used by replay.py."""
import csv
import gzip
import hashlib
import io
import json
import time
from urllib.error import URLError
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
WINDOWS = {
    "prediction": ("2050-03-11", "2050-03-14", [7, 10, 399, 701, 702, 799]),
    "historical": ("2007-05-04 13:00", "2007-05-04 22:00", [10, 399, 702, 704]),
}


def fetch(job):
    label, start, end, body = job
    params = dict(format="json", COMMAND=f"'{body}'", CENTER="'500@0'",
                  MAKE_EPHEM="'YES'", EPHEM_TYPE="'VECTORS'", START_TIME=f"'{start}'",
                  STOP_TIME=f"'{end}'", STEP_SIZE="'10 m'", OUT_UNITS="'KM-S'",
                  VEC_TABLE="'2'", CSV_FORMAT="'YES'", REF_PLANE="'FRAME'",
                  VEC_CORR="'NONE'", OBJ_DATA="'YES'")
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(params)
    for attempt in range(4):
        try:
            with urlopen(url, timeout=120) as reply:
                result = json.load(reply)["result"]
            break
        except URLError:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    header, data = result.split("$$SOE")
    data = data.split("$$EOE")[0]
    source = next(line.strip() for line in header.splitlines() if "Target body name:" in line)
    return label, body, data, {"url": url, "target": source, "rows": len(data.strip().splitlines())}


def main():
    jobs = [(label, start, end, body) for label, (start, end, bodies) in WINDOWS.items()
            for body in bodies]
    manifest = {"frame": "ICRF", "units": "TDB Julian days, km, km/s", "requests": {}}
    buffers = {label: io.StringIO() for label in WINDOWS}
    writers = {label: csv.writer(buffer, lineterminator="\n") for label, buffer in buffers.items()}
    for wr in writers.values():
        wr.writerow(["body", "jd_tdb", "x", "y", "z", "vx", "vy", "vz"])
    with ThreadPoolExecutor(max_workers=5) as pool:
        for label, body, raw, meta in pool.map(fetch, jobs):
            for row in csv.reader(raw.strip().splitlines()):
                writers[label].writerow([body, row[0].strip(), *[x.strip() for x in row[2:8]]])
            manifest["requests"][f"{label}/{body}"] = meta
    for label, buffer in buffers.items():
        encoded = buffer.getvalue().encode()
        (ROOT / "inputs" / f"{label}.csv.gz").write_bytes(gzip.compress(encoded, mtime=0))
        manifest[label + "_csv_sha256"] = hashlib.sha256(encoded).hexdigest()
    (ROOT / "horizons.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
