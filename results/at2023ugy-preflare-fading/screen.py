"""Rank pre-discovery infrared excursions in the published TDECat sample."""
import argparse
import csv
import hashlib
import io
import itertools
import math
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path
from statistics import median
from measure import EPOCH, write_csv

COMMIT = "bca013a266cd5d2a120aa4c6b34a2dc54bee1e5e"
BASE = f"https://raw.githubusercontent.com/dlangis/TDECat/{COMMIT}/"
FLARES = set("""AT2023ugy AT2023cvb AT2022upj AT2022fpx AT2022agi AT2022dyt AT2022aee AT2021uqv
AT2020afhd AT2020ksf AT2020mot AT2020nov AT2020pj AT2019qiz AT2019dsg AT2019azh AT2017gge ASASSN-14li""".split())

def get_bytes(path, remote, digest):
    if path.exists():
        data = path.read_bytes()
    else:
        with urllib.request.urlopen(BASE + urllib.parse.quote(remote), timeout=120) as response:
            data = response.read()
    actual = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    if actual != digest:
        raise ValueError(f"Wrong Git blob digest for {path.name}")
    if not path.exists():
        path.write_bytes(data)
    return data

def discovery_dates(data):
    dates = {}
    for row in csv.DictReader(io.StringIO(data.decode())):
        discovery = (date.fromisoformat(row["Discovery date (UT)"][:10]) - EPOCH).days
        for field in ["AT name", "ZTF name", "Alternative name", "Discovery internal name"]:
            if row[field]:
                dates[row[field].replace(" ", "")] = discovery
    dates["SwiftJ2058+0516"] = dates["Swift_J2058+08"]
    return dates

def band_visits(exposures):
    groups = []
    for t, mag in sorted(exposures):
        if not groups or t - groups[-1][-1][0] > 30:
            groups.append([])
        groups[-1].append((t, mag))
    return [(median(t for t, m in g), median(m for t, m in g))
            for g in groups if len(g) >= 5]

def evaluate(series, discovery):
    bands = defaultdict(list)
    for row in csv.DictReader(io.StringIO(series.decode()), delimiter=";"):
        if row["Facility"] != "NEOWISE" or row["Filter"] not in ("WISE(W1)", "WISE(W2)"):
            continue
        try:
            t, m, err = (float(row[k]) for k in ["MJD", "Magnitude", "Error"])
        except ValueError:
            continue
        if m > 0 and 0 < err < 0.5:
            bands[row["Filter"]].append((t, m))
    first = band_visits(bands["WISE(W1)"])
    second = band_visits(bands["WISE(W2)"])
    pairs = [(t, m, n) for t, m in first for u, n in second
             if abs(t-u) <= 5 and t <= discovery - 180]
    if len(pairs) < 8 or pairs[-1][0] - pairs[0][0] < 1000:
        return len(pairs), None
    groups = defaultdict(list)
    for idx, (t, m, n) in enumerate(pairs):
        direction = math.cos(2 * math.pi * (t - pairs[0][0]) / 365.25) > 0
        groups[direction].append(idx)
    if len(groups) < 2 or min(map(len, groups.values())) < 2:
        return len(pairs), None
    adjusted = [[0., 0.] for _ in pairs]
    for ids in groups.values():
        for band in (0, 1):
            offset = median(pairs[i][1+band] for i in ids)
            for i in ids:
                adjusted[i][band] = pairs[i][1+band] - offset
    windows = [[(a+b)/2 for a,b in zip(x,y)] for x,y in zip(adjusted,adjusted[1:])]
    score = max(min(x[0]-y[0], x[1]-y[1]) for x,y in itertools.permutations(windows, 2))
    return len(pairs), score

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.cache.mkdir(exist_ok=True, parents=True)
    args.out.mkdir(exist_ok=True, parents=True)
    data = get_bytes(args.cache/"catalog.csv", "Section 1: TDECat/TDE_catalogue_all.csv",
                     "6ff1376b7e5a0118a67b92cc64943c0173db119e")
    dates = discovery_dates(data)
    manifest = Path(__file__).with_name("photometry_manifest.csv")
    with manifest.open() as handle:
        objects = list(csv.DictReader(handle))
    def process(obj):
        name = obj["name"]
        path = f"Section 2: Photometry/Optical and Infrared/target_{name}_photometry.csv"
        raw = get_bytes(args.cache/(name+".csv"), path, obj["git_blob_sha1"])
        n, score = evaluate(raw, dates[name.replace(" ", "")])
        return dict(name=name, paired_visits=n, eligible=score is not None,
                    excursion_mag=score, published_ir_flare=name in FLARES)
    with ThreadPoolExecutor(max_workers=4) as pool:
        result = list(pool.map(process, objects))
    result.sort(key=lambda r: (r["excursion_mag"] is None, -(r["excursion_mag"] or 0), r["name"]))
    write_csv(args.out/"screen.csv", result)

if __name__ == "__main__":
    main()
