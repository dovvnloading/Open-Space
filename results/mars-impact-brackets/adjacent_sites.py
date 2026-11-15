"""Measure two scar centers in projected HiRISE imagery and show their CTX histories.

Usage: python adjacent_sites.py DIRECTORY OUTPUT.png
Downloads the three pinned CTX strips. Place HiRISE images and labels in DIRECTORY.
"""
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import gaussian_filter

from compare_ctx import ROOT, grayscale, raw_strip


def center(directory, product, seed, sigma):
    image = np.asarray(Image.open(directory / (product + "_RED.browse.jpg")))
    text = (directory / (product + "_RED.LBL")).read_text()

    def key(name):
        return float(re.search(r"^\s*" + name + r"\s*=\s*([-\d.]+)", text, re.M)[1])

    # A fixed local window isolates the manually identified diffuse blast zone.
    x, y = seed
    cut = image[y-150:y+150, x-150:x+150].astype(float)
    smoothed = gaussian_filter(cut, sigma)[100:200, 100:200]
    v, u = np.unravel_index(np.argmin(smoothed), smoothed.shape)
    x, y = x + u - 50, y + v - 50
    res = key("MAP_RESOLUTION")
    lon = key("CENTER_LONGITUDE") + (
        x * key("LINE_SAMPLES") / image.shape[1] - key("SAMPLE_PROJECTION_OFFSET") - 1
    ) / (res * math.cos(math.radians(key("CENTER_LATITUDE"))))
    lat = (key("LINE_PROJECTION_OFFSET") + 1 - y * key("LINES") / image.shape[0]) / res
    return float(lat), float(lon)


def separation(a, b):
    lat, lon, other_lat, other_lon = map(math.radians, (*a, *b))
    hav = math.sin((lat-other_lat)/2)**2 + math.cos(lat)*math.cos(other_lat)*math.sin((lon-other_lon)/2)**2
    return 2 * 3389.5 * math.asin(math.sqrt(hav))


def main(directory, destination):
    for sigma in [8, 12, 16, 20]:
        east = center(directory, "ESP_034503_1780", (900, 1330), sigma)
        west = center(directory, "ESP_034859_1780", (972, 1335), sigma)
        print(f"sigma={sigma}, east={east}, west={west}, separation={separation(east, west):.3f} km")

    # Fixed cutouts centered by matching local terrain. Each row is a single site.
    # The Dec 2013 eastern and western marks are 418 native CTX samples apart.
    points = [[(3683,506), (3564,487), (1256,463)], [(3273,402), (3177,410), (851,359)]]
    dates = ["18 Feb 2012", "12 Jun 2012", "6 Dec 2013"]
    sheet = Image.new("RGB", (750, 484), "white")
    draw = ImageDraw.Draw(sheet)
    with (ROOT / "adjacent-strips.csv").open() as f:
        for col, spec in enumerate(csv.DictReader(f)):
            arr = raw_strip(spec)
            for row, name in enumerate(["Eastern scar", "Western scar"]):
                x, y = points[row][col]
                patch = grayscale(arr[y-85:y+85, x-100:x+100])
                top = 242*row
                draw.text((250*col+6, top+4), f"{name}   {dates[col]}", fill="black")
                sheet.paste(patch.resize((250,212), Image.Resampling.BILINEAR), (250*col, top+24))
    sheet.save(destination)


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
