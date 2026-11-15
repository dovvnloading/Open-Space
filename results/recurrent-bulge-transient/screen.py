"""Rebuild the cross-year microlensing alert screen from the OGLE observing-season pages."""
import argparse
import csv
import hashlib
from pathlib import Path
import re
from urllib.request import urlopen

from lxml import html
import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
SNAPSHOT = HERE / "alert-seasons.csv"


def degrees(value):
    h, m, s = map(float, value.lstrip("+-").split(":"))
    return (h + m/60 + s/3600) * (-1 if value.startswith("-") else 1)


def entries(page):
    for tr in html.fromstring(page).xpath("//tr"):
        cells = [" ".join(td.text_content().split()) for td in tr.xpath("./td")]
        name = next((x for x in cells if re.fullmatch(r"\d{4}-(?:BLG|GD|DG)-\d+", x, re.I)), "")
        if not name:
            continue
        idx = next((j for j, x in enumerate(cells) if re.fullmatch(r"\d{2}:\d{2}:\d{2}(?:\.\d+)?", x)), None)
        if idx is not None:
            yield name, (15 * degrees(cells[idx]), degrees(cells[idx+1]))


def run(cache, output):
    cache.mkdir(exist_ok=True, parents=True)
    records = {}
    with SNAPSHOT.open() as handle:
        for row in csv.DictReader(handle):
            p = cache / (row["year"] + ".html")
            if not p.exists():
                p.write_bytes(urlopen(row["url"], timeout=50).read())
            page = p.read_bytes()
            if hashlib.sha256(page).hexdigest() != row["sha256"]:
                print(f'Live page differs from snapshot for {row["year"]}')
            records.update(entries(page))
    names = list(records)
    ra, dec = np.radians(np.array(list(records.values())).T)
    vec = np.column_stack((np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)))
    pairs = cKDTree(vec).query_pairs(2*np.sin(np.radians(2/3600)/2))
    candidates = sorted((np.degrees(2*np.arcsin(np.linalg.norm(vec[i]-vec[j])/2))*3600,
                         names[i], names[j]) for i, j in pairs if names[i][:4] != names[j][:4])
    with output.open("w", newline="") as dest:
        w = csv.writer(dest, lineterminator="\n")
        w.writerow(["separation_arcsec", "first_alert", "second_alert"])
        for sep, first, second in candidates:
            w.writerow([f"{sep:.6f}", first, second])
    print(f"{len(records)} positions, {len(candidates)} cross-year pairs within 2 arcsec")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("alerts-cache"))
    parser.add_argument("--out", type=Path, default=Path("pairs.csv"))
    options = parser.parse_args()
    run(options.cache, options.out)
