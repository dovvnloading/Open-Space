"""Measure pre-outburst variability of the AT 2023ugy host against the local field."""
import argparse
import csv
import hashlib
import json
import math
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import date, timedelta
from pathlib import Path
from statistics import median

TARGET = (325.18341667, 21.01626944)
EPOCH = date(1858, 11, 17)
EARLY = (2014, 2015)
LATE = (2021, 2022)
IRSA = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
GAIA = "https://gsaweb.ast.cam.ac.uk/alerts/alert/Gaia23ctp/lightcurve.csv"
FIELDS = {
    "allwise.tbl": ("allwise_p3as_psd", "ra,dec,designation,w1mpro,w2mpro"),
    "neowise.tbl": ("neowiser_p1bs_psd", "ra,dec,mjd,w1mpro,w1sigmpro,w2mpro,w2sigmpro,cc_flags,qual_frame,saa_sep,moon_masked,ph_qual,scan_id,frame_num,source_id"),
}

def sources():
    urls = {"gaia.csv": GAIA}
    for filename, (catalog, columns) in FIELDS.items():
        params = dict(catalog=catalog, spatial="Cone", objstr="325.18341667 21.01626944",
                      radius=96, radunits="arcsec", outfmt=1, selcols=columns)
        urls[filename] = IRSA + "?" + urllib.parse.urlencode(params)
    return urls

def ipac(path):
    columns = None
    for line in path.read_text().splitlines():
        if line.startswith("|") and columns is None:
            columns = [field.strip() for field in line.strip("|").split("|")]
        elif columns and line.startswith(" "):
            fields = line.split()
            if len(fields) != len(columns):
                raise ValueError(f"Unexpected IPAC row in {path.name}")
            yield dict(zip(columns, fields))

def separation(ra, dec, center):
    # Tangent plane on a field smaller than two arcminutes.
    return 3600 * math.hypot((ra-center[0]) * math.cos(math.radians(center[1])), dec-center[1])

def visits(rows):
    groups = []
    for row in sorted(rows, key=lambda item: float(item["mjd"])):
        mjd = float(row["mjd"])
        if not groups or mjd - float(groups[-1][-1]["mjd"]) > 30:
            groups.append([])
        groups[-1].append(row)
    keys = ("mjd", "w1mpro", "w2mpro", "ra", "dec")
    return [dict(n=len(group), **{k: median(float(r[k]) for r in group) for k in keys})
            for group in groups if len(group) >= 5]

def in_window(mjd, window):
    year = (EPOCH + timedelta(days=mjd)).year
    return window[0] <= year <= window[1]

def field_measurements(stars, detections):
    by_name = defaultdict(list)
    positions = {s["designation"]: (float(s["ra"]), float(s["dec"])) for s in stars}
    for row in detections:
        if row["cc_flags"][:2] != "00" or row["moon_masked"][:2] != "00":
            continue
        if int(row["qual_frame"]) != 10 or float(row["saa_sep"]) <= 5:
            continue
        if any(flag not in "AB" for flag in row["ph_qual"][:2]):
            continue
        matches = [name for name, pos in positions.items()
                   if separation(float(row["ra"]), float(row["dec"]), pos) <= 2.5]
        if len(matches) == 1:
            by_name[matches[0]].append(row)
    field = []
    target = min(stars, key=lambda s: float(s["dist"]))["designation"]
    target_visits = []
    for star in stars:
        name = star["designation"]
        series = visits(by_name[name])
        if name == target:
            target_visits = [dict(**v, offset_arcsec=separation(v["ra"], v["dec"], positions[name]))
                             for v in series]
        early, late = [[v for v in series if in_window(v["mjd"], win)] for win in (EARLY, LATE)]
        if min(len(early), len(late)) < 2:
            continue
        row = dict(source=name, offset_arcsec=float(star["dist"]), target=name == target,
                   n_visits=len(series), early_visits=len(early), late_visits=len(late),
                   early_exposures=sum(v["n"] for v in early), late_exposures=sum(v["n"] for v in late))
        for band in ["w1", "w2"]:
            for label, data in [("early", early), ("late", late)]:
                row[f"{band}_{label}"] = median(v[f"{band}mpro"] for v in data)
            row[f"{band}_fading"] = row[f"{band}_late"] - row[f"{band}_early"]
        field.append(row)
    return field, target_visits

def gaia_measurements(path):
    nights = defaultdict(list)
    for row in csv.reader(path.read_text().splitlines()[2:]):
        try:
            mjd = float(row[1]) - 2400000.5
            mag = float(row[2])
        except (ValueError, IndexError):
            continue
        nights[int(mjd)].append(mag)
    out = {}
    for label, window in [("early", EARLY), ("late", LATE)]:
        included = [mags for day, mags in nights.items() if in_window(day, window)]
        out[label] = dict(n_nights=len(included), n_transits=sum(map(len, included)),
                          median_G=median(median(g) for g in included))
    return out

def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        for row in rows:
            writer.writerow({k: round(v, 8) if isinstance(v, float) else v for k, v in row.items()})

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    receipt = {}
    for name, url in sources().items():
        path = args.cache / name
        if not path.exists():
            with urllib.request.urlopen(url, timeout=180) as response:
                raw = response.read()
            if not raw.startswith((b"\\", b"Gaia23ctp")):
                raise ValueError(f"Archive did not return the expected product for {name}")
            path.write_bytes(raw)
        raw = path.read_bytes()
        receipt[name] = {"url": url, "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}
    stars = list(ipac(args.cache / "allwise.tbl"))
    detections = list(ipac(args.cache / "neowise.tbl"))
    field, series = field_measurements(stars, detections)
    write_csv(args.out / "field.csv", field)
    write_csv(args.out / "visits.csv", series)
    result = {"position_deg": TARGET, "allwise_sources": len(stars), "neowise_detections": len(detections),
              "gaia_G": gaia_measurements(args.cache / "gaia.csv"), "inputs": receipt}
    (args.out / "measurement.json").write_text(json.dumps(result, indent=2) + "\n")

if __name__ == "__main__":
    main()
