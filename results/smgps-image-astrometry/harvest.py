"""Fetch one block of original FITS metadata per file and the image index."""

import csv
import hashlib
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

INDEX = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2BA/693/A247/list&-out=Name,GLON,GLAT,Nx,Ny,FileName&-out.max=unlimited"
BASE = "https://archive-gw-1.kat.ac.za/public/repository/10.48479/0n8c-5q84/data/SMGPS_fits/"
HERE = Path(__file__).resolve().parent


def get(url, **headers):
    for attempt in range(4):
        try:
            with urlopen(Request(url, headers=headers), timeout=40) as resp:
                return resp.read(11520 if "Range" in headers else 8000000)
        except OSError:
            if attempt == 3:
                raise
            time.sleep(2**attempt)


def filename_url(filename):
    return BASE + quote(filename, safe="")


def read_index():
    data = get(INDEX)
    rows = []
    for line in data.decode().splitlines():
        words = [w.strip() for w in line.split("\t")]
        if len(words) == 6 and words[-1].endswith(".fits"):
            name, lon, lat, nx, ny, filename = words
            rows.append(
                {
                    "name": name,
                    "filename": filename,
                    "listed_l": float(lon),
                    "listed_b": float(lat),
                    "nx": int(nx),
                    "ny": int(ny),
                }
            )
    if not rows or len({r["filename"] for r in rows}) != len(rows):
        raise ValueError("Image index is empty or has duplicate filenames")
    return rows


def add_header(row):
    filename = row["filename"]
    cache = HERE / "source-cache" / filename
    if cache.exists():
        data = cache.read_bytes()
    else:
        data = get(filename_url(filename), Range="bytes=0-11519")
        if data[:8] != b"SIMPLE  ":
            raise ValueError(f"Unexpected FITS response for {filename}")
        cache.write_bytes(data)
    header = {}
    for i in range(0, len(data), 80):
        card = data[i : i + 80].decode("ascii")
        keyword = card[:8].strip()
        if keyword == "END":
            break
        if card[8:10] == "= ":
            header[keyword] = card[10:].split("/", 1)[0].strip().strip("'").strip()
    else:
        raise ValueError(f"Truncated FITS header for {filename}")
    fields = [
        "NAXIS1",
        "NAXIS2",
        "CRPIX1",
        "CRPIX2",
        "CRVAL1",
        "CRVAL2",
        "CDELT1",
        "CDELT2",
        "CTYPE1",
        "CTYPE2",
        "CROTA1",
        "CROTA2",
    ]
    for key in fields:
        row[key] = header[key]
    extra = [
        k
        for k in header
        if k.startswith(
            ("PC1_", "PC2_", "CD1_", "CD2_", "PV1_", "PV2_", "LONPOLE", "LATPOLE")
        )
    ]
    if extra:
        raise ValueError(f"Additional projection keywords in {filename}: {extra}")
    row["header_sha256"] = hashlib.sha256(data[: i + 80]).hexdigest()
    return row


def main():
    (HERE / "source-cache").mkdir(exist_ok=True)
    rows = read_index()
    with ThreadPoolExecutor(max_workers=7) as pool:
        rows = list(pool.map(add_header, rows))
    with (HERE / "headers.csv").open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} native headers retrieved")


if __name__ == "__main__":
    main()
