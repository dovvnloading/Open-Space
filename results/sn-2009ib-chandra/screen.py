"""Optional survey-wide replay. Downloads the full CSC and optical list if needed."""
import argparse
import csv
import gzip
from pathlib import Path
import re
from urllib.request import urlretrieve

import numpy as np
from scipy.spatial import cKDTree
from lxml import html


def acquire(path, url):
    if not path.is_file():
        urlretrieve(url, path)
    return path


def optical(cache):
    positions = {}
    for year in range(2009, 2027):
        folder = f"sn{year}" if year < 2026 else "snimages"
        path = acquire(cache / f"sn{year}.html",
                       f"https://www.rochesterastronomy.org/{folder}/snredshift.html")
        table = max(html.parse(path).xpath("//table"), key=lambda t: len(t.xpath("./tr")))
        for row in table.xpath("./tr")[1:]:
            cells = [c.text_content().strip() for c in row.xpath("td")]
            if len(cells) != 16 or not re.match(r"I[aIbVc]|SLSN|TDE", cells[1]):
                continue
            def angle(s):
                h, m, sec = map(float, s.lstrip("+-").split(":"))
                return (h + m/60 + sec/3600) * (-1 if s.startswith("-") else 1)
            positions.setdefault(cells[0], [15*angle(cells[3]), angle(cells[4])])
    return np.deg2rad(np.array(list(positions.values())))


def csc(cache):
    path = acquire(cache / "csc.gz", "https://heasarc.gsfc.nasa.gov/FTP/heasarc/"
                   "dbase/tdat_files/heasarc_csc.tdat.gz")
    points = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("<DATA>"):
                break
        for row in csv.reader(fh, delimiter="|"):
            if row[0] == "<END>":
                break
            points.append([float(row[1]), float(row[2])])
    return np.deg2rad(np.array(points))


def basis(coords):
    a, d = coords.T
    c, s = np.cos(a), np.sin(a)
    radial = np.column_stack((c*np.cos(d), s*np.cos(d), np.sin(d)))
    east = np.column_stack((-s, c, np.zeros_like(s)))
    north = np.column_stack((-c*np.sin(d), -s*np.sin(d), np.cos(d)))
    return radial, east, north


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    coords = optical(args.cache)
    xray = csc(args.cache)
    sky, east, north = basis(coords)
    tree = cKDTree(basis(xray)[0])
    measured = np.deg2rad(0.34190304 / 3600)
    radius = 2*np.sin(measured/2)
    print("Optical positions, X-ray sources:", len(coords), len(xray))
    original = tree.query(sky)[0]
    print("Matches within measured offset, 1.4 arcsec:",
          np.sum(original <= radius), np.sum(original <= 2*np.sin(np.deg2rad(1.4/3600)/2)))
    for shift in [10, 20, 30]:
        h = np.deg2rad(shift/3600)
        # Exclude points too far away to reach an X-ray source under any offset.
        near = original < 2*np.sin((h+measured)/2) + 1e-12
        counts = []
        for theta in np.linspace(0, 2*np.pi, 360, endpoint=False):
            moved = (np.cos(h)*sky[near] + np.sin(h)*
                     (east[near]*np.cos(theta) + north[near]*np.sin(theta)))
            counts.append(int(np.sum(tree.query(moved)[0] <= radius)))
        print("Shift arcsec, directions, total matches, fraction with >=1:",
              shift, len(counts), sum(counts), sum(c >= 1 for c in counts)/len(counts))


if __name__ == "__main__":
    main()
