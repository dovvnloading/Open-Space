"""Replay the original map and offset-position search on the archived sample."""
import argparse
import csv
import io
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from measure import CUTOUT, ROOT


def replay(row):
    query = urlencode({"ID": row["uri"], "CIRCLE": f"{row['ra']} {row['dec']} 0.0125"})
    with urlopen(CUTOUT + "?" + query, timeout=90) as response:
        raw = response.read()
    with fits.open(io.BytesIO(raw)) as hdus:
        header = hdus[0].header
        data = hdus[0].data.squeeze().astype(float) * 1000
        wcs = WCS(header).celestial
    x, y = wcs.world_to_pixel_values(float(row["ra"]), float(row["dec"]))
    yy, xx = np.indices(data.shape)
    dx, dy = xx - x, yy - y
    radius = np.hypot(dx, dy)
    sky = data[(radius > 10) & (radius < 40) & np.isfinite(data)]
    bg = np.median(sky)
    noise = 1.482602218 * np.median(np.abs(sky - bg))
    center = np.where(radius < 3.5, data, np.nan)
    py, px = np.unravel_index(np.nanargmax(center), center.shape)
    controls = []
    for a in range(-24, 25, 8):
        for b in range(-24, 25, 8):
            if 64 <= a*a + b*b <= 28**2:
                mask = (dx-a)**2 + (dy-b)**2 < 3.5**2
                values = data[mask & np.isfinite(data)]
                if values.size:
                    controls.append((np.max(values) - bg) / noise)
    scale = abs(wcs.pixel_scale_matrix[0,0]) * 3600
    return {**row, "utc": header['DATE-OBS'][:19],
            "peak_mJy_beam": round(data[py,px],3),
            "central_pixel_mJy_beam": round(data[round(float(y)),round(float(x))],3),
            "image_rms_mJy_beam": round(noise,3),
            "peak_snr": round((data[py,px]-bg)/noise,2),
            "peak_offset_as": round(radius[py,px]*scale,3),
            "n_controls": len(controls), "control_max_snr":round(max(controls),2),
            "controls_ge10":sum(v>=10 for v in controls),
            "controls_ge18":sum(v>=18 for v in controls),
            "controls_ge22":sum(v>=22 for v in controls)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--name",help="Replay only one object")
    args=parser.parse_args()
    if args.output.resolve() == (ROOT/"screen.csv").resolve():
        parser.error("Choose a separate path for the replay output")
    with (ROOT/"screen.csv").open() as source, args.output.open('w') as target:
        rows=csv.DictReader(source)
        out=csv.DictWriter(target,rows.fieldnames,lineterminator="\n")
        out.writeheader()
        for row in rows:
            if args.name is None or row['name']==args.name:
                out.writerow(replay(row))
                target.flush()


if __name__ == "__main__":
    main()
