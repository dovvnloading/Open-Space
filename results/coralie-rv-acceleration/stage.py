"""Freeze a CORALIE14 observing sample and retrieve exact FITS header ranges."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timedelta, timezone
from hashlib import sha256
from io import BytesIO
import itertools
import json
from pathlib import Path
import re
import threading
import time

from astropy.io import fits
from astropy.table import Table
import numpy as np
import requests


EPOCH = datetime(1858, 11, 17)
WINDOWS = {
    "training": (datetime(2015, 7, 1), datetime(2020, 1, 1)),
    "holdout": (datetime(2020, 1, 1), datetime(2024, 10, 15)),
}
DPID = re.compile(r"ADP\.\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}")
ARCHIVE = "https://dataportal.eso.org/dataPortal/file/"
TLS = threading.local()
COLUMNS = (
    "target", "split", "year", "night", "id", "mjd", "ra", "dec", "snr",
)
HEADER_COLUMNS = (
    "object", "obs_mjd", "obs_ra", "obs_dec", "rvmps", "errmps", "mask", "qc",
    "bismps", "fwhmmps", "ref", "fp", "provenance", "pipeline",
)
RECEIPT_COLUMNS = (*COLUMNS, *HEADER_COLUMNS, "bytes", "sha256", "http_range", "outcome")


def write_csv(path, columns, records):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def chunks_by_night(rows):
    by_night = {}
    for row in rows:
        current = by_night.get(row["night"])
        if current is None or (row["snr"], row["id"]) > (current["snr"], current["id"]):
            by_night[row["night"]] = row
    return sorted(by_night.values(), key=lambda row: row["mjd"])


def spaced_subsample(rows):
    if len(rows) <= 18:
        return rows
    selected = set()
    for _, group in itertools.groupby(range(len(rows)), key=lambda i: rows[i]["year"]):
        members = list(group)
        if len(members) <= 3:
            selected.update(members)
        else:
            selected.update((members[0], members[len(members) // 2], members[-1]))
    while len(selected) < 18:
        remaining = set(range(len(rows))) - selected
        chosen = max(remaining, key=lambda i: (
            min(abs(rows[i]["mjd"] - rows[j]["mjd"]) for j in selected), -i,
        ))
        selected.add(chosen)
    return [rows[i] for i in sorted(selected)]


def as_string(value):
    return value.decode().strip() if isinstance(value, bytes) else str(value).strip()


def prepare(metadata, paper_dir, destination):
    tap = json.loads((paper_dir / "tap-receipt.json").read_text(encoding="utf-8"))
    if sha256(metadata.read_bytes()).hexdigest() != tap["sha256"]:
        raise ValueError("Metadata bytes differ from the frozen ESO receipt")
    tab = Table.read(metadata, format="votable")
    windows = {key: tuple((item - EPOCH).days for item in value) for key, value in WINDOWS.items()}
    targets = {}
    for row in tab:
        mjd = float(row["t_min"])
        phase = next((key for key, (a, b) in windows.items() if a <= mjd < b), None)
        if phase is None:
            continue
        name = as_string(row["target_name"]).upper()
        ident = as_string(row["dp_id"])
        if not DPID.fullmatch(ident):
            raise ValueError(f"Unexpected dataset ID {ident!r}")
        snr = float(row["snr"]) if not np.ma.is_masked(row["snr"]) else -1
        entry = dict(target=name, split=phase, year=(EPOCH + timedelta(days=mjd)).year,
                     night=int(mjd - 0.5), id=ident, mjd=mjd,
                     ra=float(row["s_ra"]), dec=float(row["s_dec"]), snr=snr)
        targets.setdefault(name, {key: [] for key in WINDOWS})[phase].append(entry)

    selected = []
    coverage = []
    for name, grouped in sorted(targets.items()):
        nights = {key: chunks_by_night(value) for key, value in grouped.items()}
        counts = {key: len(value) for key, value in nights.items()}
        years = {key: len({row["year"] for row in value}) for key, value in nights.items()}
        if all(counts[key] >= 12 and years[key] >= 4 for key in WINDOWS):
            for key, items in nights.items():
                selected.extend(spaced_subsample(items))
            coverage.append({"target": name, **{f"{key}_nights": value for key, value in counts.items()},
                             **{f"{key}_years": value for key, value in years.items()}})
    if len({row["id"] for row in selected}) != len(selected):
        raise ValueError("Duplicate dataset IDs in selected sample")
    write_csv(destination, COLUMNS, selected)
    write_csv(destination.with_name("coverage.csv"), tuple(coverage[0]), coverage)
    print(f"{len(tab)} archive spectra, {len(targets)} names in date windows, "
          f"{len(coverage)} eligible names, {len(selected)} fixed spectra", flush=True)


def header_fields(header, item):
    if header["ARCFILE"] != item["id"] + ".fits" or header["INSTRUME"] != "CORALIE":
        raise ValueError(f"Wrong instrument/product in {item['id']}")
    if abs(float(header["MJD-OBS"]) - float(item["mjd"])) > 2 / 86400:
        raise ValueError(f"FITS MJD differs from TAP in {item['id']}")

    def to_ms(key):
        value = header.get("ESO QC CCF " + key)
        return None if value is None else 1000 * float(value)

    return dict(object=header.get("OBJECT"), obs_mjd=header["MJD-OBS"], obs_ra=header.get("RA"),
                obs_dec=header.get("DEC"), rvmps=to_ms("RV"), errmps=to_ms("RV ERROR"),
                bismps=to_ms("BIS SPAN"), fwhmmps=to_ms("FWHM"),
                mask=header.get("ESO QC CCF MASK"), qc=header.get("ESO QC SCIRED CHECK"),
                ref=header.get("ESO PRO REC1 RAW1 CATG"), fp=header.get("ESO QC FIB B FP"),
                provenance=header.get("ESO PRO REC1 RAW1 NAME"), pipeline=header.get("PROCSOFT"))


def session():
    if not hasattr(TLS, "client"):
        TLS.client = requests.Session()
        TLS.client.headers["User-Agent"] = "Open-Space-CORALIE/1.0 scientific-archive-replay"
    return TLS.client


def fetch_range(ident, size):
    with session().get(ARCHIVE + ident, headers={"Range": f"bytes=0-{size - 1}"},
                       stream=True, timeout=(12, 30)) as response:
        response.raise_for_status()
        content_range = response.headers.get("Content-Range", "")
        if response.status_code != 206 or not re.fullmatch(fr"bytes 0-{size - 1}/\d+", content_range):
            raise ValueError(f"Unexpected HTTP response: {response.status_code} {content_range}")
        block = response.content
        if len(block) != size:
            raise ValueError(f"Expected {size} response bytes")
    return block, content_range


def retrieve(item, data, offline):
    ident = item["id"]
    if not DPID.fullmatch(ident):
        raise ValueError(f"Unexpected dataset ID {ident!r}")
    path = data / item["split"] / (ident.replace(":", "_") + ".header")
    receipt = path.with_suffix(".range")
    try:
        if path.is_file() and receipt.is_file():
            block, content_range = path.read_bytes(), receipt.read_text(encoding="ascii")
            header = fits.Header.fromfile(BytesIO(block))
        elif offline:
            return dict(item, outcome="missing staged bytes")
        else:
            try:
                block, content_range = fetch_range(ident, 115200)
            except requests.RequestException:
                time.sleep(1.25)
                block, content_range = fetch_range(ident, 115200)
            try:
                header = fits.Header.fromfile(BytesIO(block))
            except OSError:
                block, content_range = fetch_range(ident, 230400)
                header = fits.Header.fromfile(BytesIO(block))
            path.write_bytes(block)
            receipt.write_text(content_range, encoding="ascii")
        return dict(item, **header_fields(header, item), bytes=len(block), sha256=sha256(block).hexdigest(),
                    http_range=content_range, outcome="ok")
    except (requests.RequestException, OSError, KeyError, ValueError) as exc:
        return dict(item, outcome=f"{type(exc).__name__}: {exc}")


def run_download(selection, data, output, workers, offline):
    with selection.open(encoding="utf-8", newline="") as stream:
        chosen = list(csv.DictReader(stream))
    phases = {row["split"] for row in chosen}
    if not phases <= {*WINDOWS, "coralie98", "coralie07"}:
        raise ValueError("Unknown observing phase in the selection")
    for phase in phases:
        (data / phase).mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = []
        for index, row in enumerate(pool.map(lambda item: retrieve(item, data, offline), chosen), start=1):
            results.append(row)
            if index % 100 == 0:
                succeeded = sum(r["outcome"] == "ok" for r in results)
                print(f"{index}/{len(chosen)} processed, {succeeded} available", flush=True)
    write_csv(output, RECEIPT_COLUMNS, results)
    print(f"{sum(r['outcome'] == 'ok' for r in results)}/{len(chosen)} complete header receipts", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    prepare_parser = commands.add_parser("freeze")
    prepare_parser.add_argument("metadata", type=Path)
    prepare_parser.add_argument("selection", type=Path)
    download_parser = commands.add_parser("stage")
    download_parser.add_argument("selection", type=Path)
    download_parser.add_argument("data", type=Path)
    download_parser.add_argument("output", type=Path)
    download_parser.add_argument("--workers", type=int, default=12)
    download_parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    if args.action == "freeze":
        prepare(args.metadata, Path(__file__).resolve().parent, args.selection)
    else:
        run_download(args.selection, args.data, args.output, args.workers, args.offline)


if __name__ == "__main__":
    main()
