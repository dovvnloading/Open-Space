"""Verify 76 published positions against the original HiRISE RED product."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import tempfile
import zipfile

import numpy as np
from openpyxl import load_workbook
from PIL import Image, ImageDraw
from scipy.ndimage import gaussian_laplace, map_coordinates
from scipy.stats import spearmanr

IMAGE_ID = "ESP_066902_1830"
ARCHIVE_SHA = "f32fd7ece38d8ebb9a8874651d14f09c450278f8e58a27858c12f88219651aed"
JP2_SHA = "93bf4eca971b67b9f88561127e7804c2cf0a8036b6049e9ce54dc48144fbe0c1"
LABEL_SHA = "44feb30aabfb4a55de07a23298692e64b3f8e82fe231ddaeedf7863e7e017e28"
WINDOW = (10200, 17700, 1050, 1100)


def verify(path, expected):
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if digest != expected:
        raise ValueError(f"Input checksum mismatch for {path}")


def number(label, key):
    return float(re.search(rf"^\s*{key}\s*=\s*([-+0-9.]+)", label, re.M).group(1))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("archive", "jp2", "label", "output"):
        p.add_argument(key, type=Path)
    args = p.parse_args()
    for path, expected in [(args.archive, ARCHIVE_SHA), (args.jp2, JP2_SHA), (args.label, LABEL_SHA)]:
        verify(path, expected)
    with zipfile.ZipFile(args.archive) as z:
        data = z.read(f"cluster_data/data_set/{IMAGE_ID}.xlsx")
    rows = list(load_workbook(io.BytesIO(data), data_only=True, read_only=True).active.values)
    diameter, lon, lat = np.array([[row[rows[0].index(key)] for key in ("Diam_km", "x_coord", "y_coord")]
                                   for row in rows[1:]], dtype=float).T
    label = args.label.read_text()
    resolution = number(label, "MAP_RESOLUTION")
    scale = number(label, "MAP_SCALE")
    x = ((lon % 360) - number(label, "CENTER_LONGITUDE")) * resolution
    x += number(label, "SAMPLE_PROJECTION_OFFSET") - WINDOW[0]
    y = number(label, "LINE_PROJECTION_OFFSET") - lat * resolution - WINDOW[1]
    radius = 500 * diameter / scale
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "cut.png"
        subprocess.run(["gdal_translate", "-q", "-of", "PNG", "-ot", "UInt16", "-srcwin",
                        *map(str, WINDOW), str(args.jp2), str(path)], check=True)
        data = np.array(Image.open(path), dtype=float)
    margin = 38 + 1.8 * radius
    if not np.all((x >= margin) & (x < WINDOW[2] - margin) &
                  (y >= margin) & (y < WINDOW[3] - margin)):
        raise ValueError("An aperture would cross the image crop.")
    angles = np.arange(32) * 2 * np.pi / 32

    def central_darkening(xc, yc):
        def ring(factor):
            dx = (factor * radius[:, None]) * np.cos(angles)
            dy = (factor * radius[:, None]) * np.sin(angles)
            return map_coordinates(data, [yc[:, None] + dy, xc[:, None] + dx], order=1).mean(axis=1)
        return ring(1.8) - .5 * (ring(.3) + ring(.7))

    observed = central_darkening(x, y)
    control = [float(central_darkening(x + offset * np.cos(angle), y + offset * np.sin(angle)).mean())
               for offset in (14, 22, 38) for angle in np.arange(48) * 2 * np.pi / 48]
    report = dict(image_id=IMAGE_ID, jp2_sha256=JP2_SHA, label_sha256=LABEL_SHA,
                  workbook_zip_sha256=ARCHIVE_SHA, gdal_source_window=WINDOW, meters_per_pixel=scale,
                  positions=len(x), central_darkening_positive=int(sum(observed > 0)),
                  central_darkening_mean_DN=float(observed.mean()),
                  central_darkening_median_DN=float(np.median(observed)),
                  translated_mean_DN_min_median_max=np.quantile(control,[0,.5,1]).tolist(),
                  translated_means_DN=control)
    # Infer dark-feature width using image intensities alone. The first row
    # is the registered solution; the following 144 are translated controls.
    shifts = [(0., 0.)] + [(offset * np.cos(angle), offset * np.sin(angle))
                            for offset in (14, 22, 38) for angle in np.arange(48) * 2*np.pi/48]
    xs = np.array([x + dx for dx, _ in shifts])
    ys = np.array([y + dy for _, dy in shifts])
    scales = np.geomspace(.6, 9., 30)
    response = []
    for scale_px in scales:
        filtered = gaussian_laplace(data, scale_px) * scale_px**2
        response.append(map_coordinates(filtered, [ys, xs], order=1))
    widths = scales[np.argmax(response, axis=0)]
    correlations = [float(spearmanr(diameter, width).statistic) for width in widths]
    report["diameter_blind_log_scale"] = dict(
        scale_grid_pixels=scales.tolist(), spearman_registered=correlations[0],
        spearman_translated_min_median_max=np.quantile(correlations[1:], [0,.5,1]).tolist(),
        translated_at_least_as_high=sum(r >= correlations[0] for r in correlations[1:]),
        at_largest_grid_scale=int(sum(widths[0] == max(scales))),
        catalog_diameter_m=(1000 * diameter).tolist(), image_scale_px=widths[0].tolist(),
        spearman_translated=correlations[1:])
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "pixels.json").write_text(json.dumps(report, indent=2) + "\n")
    lo, hi = np.percentile(data, [2,98])
    image = Image.fromarray(np.uint8(255*np.clip((data-lo)/(hi-lo),0,1))).convert("RGB")
    draw = ImageDraw.Draw(image)
    for xc, yc, r in zip(x,y,radius):
        draw.ellipse((xc-r-2,yc-r-2,xc+r+2,yc+r+2), outline=(227,72,55))
    # 25 m scale on the original 0.25 m per pixel image.
    draw.rectangle((36, 978, 163, 1028), fill=(28,28,28))
    draw.line((48, 1012, 148, 1012), fill="white", width=4)
    draw.text((50, 983), "25 m", fill="white", font_size=19)
    image.save(args.output / "field.jpg", quality=87)


if __name__ == "__main__":
    main()
