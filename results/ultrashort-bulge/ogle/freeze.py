"""Freeze OGLE-IV EWS BLG catalogs and light curves, resuming verified receipts.

Run: python results/ultrashort-bulge/ogle/freeze.py CACHE [--jobs 6].
Outputs CSV inventories beneath this script; raw files stay in CACHE.
"""
import argparse
import concurrent.futures as cf
import csv
import hashlib
import json
import pathlib
import sys
from datetime import datetime, timezone

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "ogle-dark-lenses"))
from fetch import BASE, events, get


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write_csv(name, records):
    path = HERE / name
    tmp = path.with_suffix(".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    tmp.replace(path)


def catalog(cache, year):
    directory = cache / year
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "lenses.par"
    url = f"{BASE}/{year}/lenses.par"
    if not path.exists():
        data = get(url)
        if not list(events(data)):
            raise ValueError(f"empty EWS table {url}")
        path.write_bytes(data)
    data = path.read_bytes()
    names = set(events(data))
    rows = []
    for line in data.decode().splitlines():
        cells = line.split()
        if not cells or cells[0] not in names or "-BLG-" not in cells[0]:
            continue
        rows.append(dict(event=cells[0], field=cells[1], star_no=cells[2], ra_hms=cells[3],
                         dec_dms=cells[4], t0_hjd=cells[5] if len(cells)>5 else "",
                         tau_days=cells[7] if len(cells)>7 else "",
                         umin=cells[8] if len(cells)>8 else ""))
    return dict(year=year, url=url, bytes=len(data), sha256=sha(data),
                archive_ids=len(names), bulge_ids=len(rows)), rows


def fetch_curve(cache, row):
    year, _, number = row["event"].split("-")
    path = cache / year / f"blg-{number}.dat"
    url = f"{BASE}/{year}/blg-{number}/phot.dat"
    meta = dict(event=row["event"], url=url, sha256="", bytes="", n="", min_hjd="", max_hjd="",
                n_before_30d="", n_after_30d="", n_peak_0p5d="", peak_gap_hours="",
                searched_utc="", status="")
    try:
        data = path.read_bytes() if path.exists() else get(url)
        from io import BytesIO
        t = np.loadtxt(BytesIO(data), usecols=0, ndmin=1)
        if len(t) == 0 or not np.all(np.isfinite(t)) or np.any(t < 2400000):
            raise ValueError("invalid HJD column")
        if not path.exists():
            tmp = path.with_suffix(".partial")
            tmp.write_bytes(data)
            tmp.replace(path)
        meta.update(sha256=sha(data), bytes=len(data), n=len(t), min_hjd=f"{t.min():.5f}",
                    max_hjd=f"{t.max():.5f}", status="ok")
        try:
            t0 = float(row["t0_hjd"])
        except ValueError:
            t0 = np.nan
        if t.min()-1000 < t0 < t.max()+1000:
            peak = np.sort(t[np.abs(t - t0) <= 0.5])
            meta.update(n_before_30d=np.sum(t<t0-30), n_after_30d=np.sum(t>t0+30),
                        n_peak_0p5d=len(peak),peak_gap_hours=round(float(np.median(np.diff(peak))*24),3)
                        if len(peak)>1 else "")
    except Exception as exc:
        meta["status"] = f"{type(exc).__name__}: {exc}"
    meta["searched_utc"] = utc()
    return meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=pathlib.Path)
    parser.add_argument("--jobs", type=int, default=6)
    args = parser.parse_args()
    all_rows, sources = [], []
    for year in ["2025","2026"]:
        info, rows = catalog(args.cache, year)
        all_rows.extend(rows)
        sources.append(info)
    write_csv("catalog.csv", all_rows)
    (HERE/"tables.json").write_text(json.dumps(sources, indent=2)+"\n")
    def rank(row):
        try:
            te = float(row["tau_days"])
            return te if te > 0 else 5
        except ValueError:
            return 5
    out = []
    with cf.ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(fetch_curve, args.cache, row) for row in sorted(all_rows, key=rank)]
        for i, f in enumerate(cf.as_completed(futures), start=1):
            out.append(f.result())
            if i % 40 == 0 or i <= 4:
                write_csv("photometry.csv", sorted(out, key=lambda x:x["event"]))
                print(i, len(all_rows), utc(), flush=True)
    write_csv("photometry.csv", sorted(out, key=lambda x:x["event"]))


if __name__ == "__main__":
    main()
