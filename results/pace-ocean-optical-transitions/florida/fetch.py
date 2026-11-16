"""Retrieve frozen ERDDAP index windows. Save NetCDF subsets and a request ledger."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import csv
import time
from urllib.parse import quote
from urllib.request import urlopen

BASE = "https://erddap.marine.usf.edu/erddap/griddap/pace_oc_aop_1d_fl.nc?"
ETOPO = "https://upwell.pfeg.noaa.gov/erddap/griddap/etopo180.nc?"
BANDS = (413, 442, 465, 490, 510, 530, 560, 620, 665, 678)
CACHE = Path(__file__).resolve().parent / "cache"


def download(item):
    name, url = item
    path = CACHE / name
    for attempt in range(3):
        try:
            if not path.exists():
                with urlopen(url, timeout=180) as response:
                    data = response.read()
                if not data.startswith(b"CDF"):
                    raise ValueError(f"Unexpected response: {data[:80]!r}")
                partial = path.with_suffix(".tmp")
                partial.write_bytes(data)
                partial.replace(path)
            data = path.read_bytes()
            print(name, len(data), flush=True)
            return name, len(data), hashlib.sha256(data).hexdigest(), url
        except Exception as exc:
            print(name, type(exc).__name__, str(exc), flush=True)
            if attempt == 2:
                raise
            time.sleep(2**attempt)


def main():
    CACHE.mkdir(exist_ok=True)
    jobs = []
    for start in range(0, 687, 36):
        end = min(686, start + 35)
        slab = f"[{start}:4:{end}][116:24:764][12:24:1092]"
        query = ",".join(
            v + slab for v in [*(f"Rrs_{b}" for b in BANDS), "STRAYLIGHT_MASK"]
        )
        jobs.append((f"oci_{start:03}.nc", BASE + quote(query, safe=",:")))
    jobs.append(
        (
            "etopo.nc",
            ETOPO + quote("altitude[(24):2:(30.1)][(-89):2:(-79)]", safe=",:()"),
        )
    )
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(download, jobs))
    with (CACHE / "manifest.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "bytes", "sha256", "url"])
        writer.writerows(results)


if __name__ == "__main__":
    main()
