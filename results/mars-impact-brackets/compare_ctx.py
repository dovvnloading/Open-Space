"""Render before/after CTX pixels with measured control-point registration."""
import csv
import hashlib
import sys
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parent


def raw_strip(spec):
    first, last = map(int, (spec["first_byte"], spec["last_byte"]))
    request = Request(spec["url"], headers={"Range": f"bytes={first}-{last}"})
    with urlopen(request, timeout=90) as response:
        if response.status != 206:
            raise ValueError("Server did not return a byte range")
        data = response.read(last-first+1)
    if hashlib.sha256(data).hexdigest() != spec["sha256"]:
        raise ValueError("Raw strip does not match input hash")
    return np.frombuffer(data, dtype=np.uint8).reshape(-1, 5056)


def estimate_registration(before, after):
    # Similarity transform estimated on four landmarks outside the blast zone.
    design = []
    for x, y in before:
        design.extend([[x, -y, 1, 0], [y, x, 0, 1]])
    a, b, tx, ty = np.linalg.lstsq(design, after.ravel(), rcond=None)[0]
    return np.array([[a, -b], [b, a]]), np.array([tx, ty])


def grayscale(crop):
    low, high = np.percentile(crop, [1, 99])
    return Image.fromarray(np.uint8(np.clip((crop-low)/(high-low), 0, 1)*255))


def main(destination):
    with (ROOT / "ctx-strips.csv").open() as f:
        before, after = [raw_strip(spec) for spec in csv.DictReader(f)]
    before_marks = np.array([[2272,380], [2584,264], [2644,78], [2422,378]])
    after_marks = np.array([[2924,366], [3258,248], [3316,66], [3090,374]])
    matrix, offset = estimate_registration(before_marks, after_marks)
    errors = np.linalg.norm(before_marks @ matrix.T + offset - after_marks, axis=1)
    print(f"Four-landmark RMS registration residual: {np.mean(errors**2)**0.5:.1f} CTX pixels")

    cx, cy, half = 3266, 456, 105
    y, x = np.mgrid[cy-half:cy+half, cx-half:cx+half]
    mapped = np.linalg.solve(matrix, np.array([x.ravel(), y.ravel()]) - offset[:, None])
    before_crop = map_coordinates(before.astype(float), mapped[::-1], order=1).reshape(x.shape)
    after_crop = after[cy-half:cy+half, cx-half:cx+half]
    sheet = Image.new("RGB", (630, 342), "white")
    draw = ImageDraw.Draw(sheet)
    for i, (im, date) in enumerate([(before_crop, "14 February 2015"), (after_crop, "19 February 2016")]):
        sheet.paste(grayscale(im).resize((315,315), Image.Resampling.BILINEAR), (315*i,27))
        draw.text((315*i+6,6), "CTX   " + date, fill="black")
    sheet.save(destination)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
