"""Project two HiRISE browse products onto a common 3 m grid and register them.

Usage: python compare_images.py IMAGE_DIRECTORY OUTPUT.png
Requires numpy, scipy and Pillow. Input URLs and hashes are in sources.csv.
"""
import math
import re
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import map_coordinates

PRODUCTS = ["ESP_047030_1955", "ESP_047175_1955"]
# Midpoint of the two catalog positions. Rounded scale affects both images equally.
LATITUDE = 15.5725
LONGITUDE = 234.0905
METERS_PER_DEGREE = 59200


def image_on_grid(directory, product):
    label = (directory / f"{product}_RED.LBL").read_text()
    with Image.open(directory / f"{product}_RED.browse.jpg") as image:
        pixels = np.asarray(image.convert("L"), dtype=float)

    def value(name):
        return float(re.search(r"^\s*" + name + r"\s*=\s*([-\d.]+)", label, re.M)[1])

    east, south = np.meshgrid(np.arange(-240, 241, 3), np.arange(-240, 241, 3))
    lat = LATITUDE - south / METERS_PER_DEGREE
    lon = LONGITUDE + east / (METERS_PER_DEGREE * math.cos(math.radians(LATITUDE)))
    res = value("MAP_RESOLUTION")
    sample = value("SAMPLE_PROJECTION_OFFSET") + 1 + res * (
        lon - value("CENTER_LONGITUDE")) * math.cos(math.radians(value("CENTER_LATITUDE")))
    line = value("LINE_PROJECTION_OFFSET") + 1 - res * lat
    sample *= pixels.shape[1] / value("LINE_SAMPLES")
    line *= pixels.shape[0] / value("LINES")
    return map_coordinates(pixels, [line, sample], order=1)


def register(first, second):
    # Same central 111-by-111 field and same number of pixels at every shift.
    patch = first[25:-25, 25:-25]
    best = (-2, 0, 0)
    for dy in range(-20, 21):
        for dx in range(-20, 21):
            moved = second[25+dy:-25+dy, 25+dx:-25+dx]
            coefficient = float(np.corrcoef(patch.ravel(), moved.ravel())[0, 1])
            best = max(best, (coefficient, dy, dx))
    return best


def main(directory, destination):
    images = [image_on_grid(directory, p) for p in PRODUCTS]
    coefficient, dy, dx = register(*images)
    print(f"Pearson r={coefficient:.6f}, offset east={3*dx} m, south={3*dy} m")
    aligned = [images[0][25:-25, 25:-25], images[1][25+dy:-25+dy, 25+dx:-25+dx]]
    result = Image.new("RGB", (666, 378), "white")
    draw = ImageDraw.Draw(result)
    for i, (crop, product, date) in enumerate(zip(aligned, PRODUCTS, ["8 August", "19 August"])):
        image = Image.fromarray(np.uint8(np.clip(crop, 0, 255))).resize((333, 333), Image.Resampling.BILINEAR)
        result.paste(image, (333*i, 34))
        draw.text((333*i+7, 7), f"{date} 2016   {product}", fill="black")
    # Each panel spans 333 m, so 100 pixels correspond to 100 m.
    draw.line((14, 346, 114, 346), fill="white", width=3)
    draw.text((14, 329), "100 m", fill="white")
    result.save(destination)


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
