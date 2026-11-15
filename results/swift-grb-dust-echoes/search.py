"""Replay the Swift XRT pseudo-distance search on checksum-pinned event lists."""
import argparse
import csv
import hashlib
import math
from datetime import datetime
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy import units as units
from scipy.special import xlogy

ROOT = Path(__file__).resolve().parent
ARCHIVE = "https://www.swift.ac.uk/archive/reproc"
TABLE = "https://swift.gsfc.nasa.gov/archive/grb_table/fullview/"
EDGES = np.arange(np.log(20), np.log(20000), 0.025)
CONTROL_NAMES = {"050713A", "050724", "061019", "070129", "160623A", "221009A"}


def rows(filename):
    with (ROOT / filename).open() as stream:
        return list(csv.DictReader(stream))


def write_csv(path, records):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0])
        writer.writeheader()
        writer.writerows(records)


def select_sample(html_bytes):
    from lxml import html
    selected = []
    for row in html.fromstring(html_bytes).xpath("//tr[td]"):
        cells = [cell.text_content().strip() for cell in row.xpath("td")]
        if len(cells) < 24 or not cells[2].isdigit():
            continue
        name = cells[0]
        try:
            trigger = datetime.strptime(name[:6] + " " + cells[1], "%y%m%d %H:%M:%S")
            position = SkyCoord(cells[13], cells[14], unit=(units.hourangle, units.deg))
        except ValueError:
            continue
        try:
            fluence = float(cells[7])
        except ValueError:
            fluence = 0.0
        latitude = position.galactic.b.deg
        rank = fluence / max(0.05, abs(math.sin(math.radians(latitude))))
        selected.append(dict(name=name, trigger=trigger.isoformat(), id=int(cells[2]),
                             ra=position.ra.deg, dec=position.dec.deg, lat=round(latitude, 3),
                             flu=fluence, rank=round(rank, 3)))
    selected.sort(key=lambda item: -item["rank"])
    return [r for i, r in enumerate(selected) if i < 120 or r["name"] in CONTROL_NAMES]


def event_path(cache, rec):
    filename = rec["filename"]
    path = cache / (rec["name"] + "_" + filename)
    if not path.exists():
        url = f"{ARCHIVE}/{rec['obsid']}/xrt/event/{filename}"
        with urlopen(url, timeout=90) as stream:
            data = stream.read()
        if hashlib.sha256(data).hexdigest() != rec["sha256"]:
            raise ValueError(f"Archive checksum changed for {filename}")
        path.write_bytes(data)
    elif hashlib.sha256(path.read_bytes()).hexdigest() != rec["sha256"]:
        raise ValueError(f"Bad cached file {path}")
    return path


def geometry(row, hdu):
    h = hdu.header
    i, j = (hdu.columns.names.index(v) + 1 for v in ("X", "Y"))
    w = WCS(naxis=2)
    for name, prefix in [("crpix", "TCRPX"), ("crval", "TCRVL"),
                         ("cdelt", "TCDLT"), ("ctype", "TCTYP")]:
        setattr(w.wcs, name, [h[f"{prefix}{k}"] for k in (i, j)])
    x, y = w.all_world2pix(float(row["ra"]), float(row["dec"]), 1)
    utc = "2022-10-09T13:16:59" if row["name"] == "221009A" else row["trigger"]
    t0 = (Time(utc, scale="utc").tt.mjd - h["MJDREFI"] - h["MJDREFF"]) * 86400
    return x, y, abs(w.wcs.cdelt[0]) * 3600, abs(w.wcs.cdelt[1]) * 3600, t0


def photons(path, row):
    with fits.open(path) as hd:
        events = hd["EVENTS"]
        x0, y0, sx, sy, t0 = geometry(row, events)
        data = events.data
        time = data["TIME"] - t0
        dx, dy = (data["X"] - x0) * sx, (data["Y"] - y0) * sy
        radius = np.hypot(dx, dy)
        energy = data["PI"]
        use = ((data["GRADE"] <= 12) & (radius >= 50) & (radius < 550) &
               (time > 200) & (time < 432000) & (energy >= 50) & (energy <= 700))
        return np.array([time[use], dx[use], dy[use], energy[use]])


def peak(array):
    time, dx, dy, pi = array
    soft = (pi >= 70) & (pi <= 300)
    dist = 826.6 * time[soft] / (dx[soft]**2 + dy[soft]**2)
    hist = np.histogram(np.log(dist), EDGES)[0]
    best = (0., 0., 0, 0)
    for i in range(14, len(hist) - 14):
        n = hist[i-3:i+4].sum()
        m = hist[i-13:i-4].sum() + hist[i+5:i+14].sum()
        if n <= m * 7/18:
            continue
        mu = (n + m) * 7/25
        score = np.sqrt(2 * (xlogy(n, n/mu) + xlogy(m, m/(n+m-mu))))
        best = max(best, (score, np.exp((EDGES[i]+EDGES[i+1])/2), n, m))
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("xrt-cache"))
    parser.add_argument("--name", help="restrict replay to one GRB, e.g. 180728A")
    parser.add_argument("--output", type=Path, default=Path("screen-replay.csv"))
    parser.add_argument("--refresh-sample", action="store_true", help="reselect from live NASA table")
    parser.add_argument("--save-photons", type=Path, help="save each search array for follow-up")
    args = parser.parse_args()
    if args.refresh_sample:
        with urlopen(TABLE, timeout=90) as src:
            write_csv(args.output, select_sample(src.read()))
        return
    args.cache.mkdir(exist_ok=True, parents=True)
    if args.save_photons:
        args.save_photons.mkdir(exist_ok=True, parents=True)
    manifest = rows("event-files.csv")
    records = []
    for row in rows("sample.csv"):
        name = row["name"]
        if args.name and name != args.name:
            continue
        files = [r for r in manifest if r["name"] == name]
        arrays = [photons(event_path(args.cache, r), row) for r in files]
        events = np.concatenate(arrays, axis=1) if arrays else np.empty((4, 0))
        if args.save_photons:
            np.save(args.save_photons / (name + ".npy"), events)
        score, distance, n, m = peak(events)
        records.append(dict(name=name, files=len(files), selected_events=events.shape[1],
                            soft_events=np.count_nonzero((events[3]>=70) & (events[3]<=300)),
                            score=round(score, 5), distance_pc=round(distance, 3),
                            on=int(n), side=int(m)))
    if not records:
        parser.error("No matching target")
    write_csv(args.output, records)


if __name__ == "__main__":
    main()
