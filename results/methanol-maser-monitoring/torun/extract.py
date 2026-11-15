"""Extract vector marker centers from Szymczak et al. (2018), supplementary PDF."""
import argparse
import hashlib
from pathlib import Path
from urllib.request import urlopen

import fitz
import numpy as np
import pandas as pd

SOURCE = "https://repozytorium.umk.pl/bitstreams/fad49269-d41d-4bab-8601-dfd7e6581847/download"
SHA256 = "3010a9ca4c8965dae33c28324ccd6c5bb7617d151f11d11896f43658ed0573ad"
HERE = Path(__file__).resolve().parent


def extract(pdf):
    if hashlib.sha256(pdf).hexdigest() != SHA256:
        raise ValueError("Unexpected PDF hash. Check the source and inspect the figure.")
    with fitz.open(stream=pdf, filetype="pdf") as doc:
        page = doc[7]
        if "G28.011−0.426" not in page.get_text():
            raise ValueError("Target panel missing")
        draws = page.get_drawings()
    # The seven x paths contain minor and major ticks. Their third segment
    # starts on a labeled major tick, 55200 through 56400 in steps of 200 days.
    major_x = [draws[i]["items"][2][1].x for i in range(1144, 1151)]
    slope, offset = np.polyfit(np.arange(55200, 56401, 200), major_x, 1)
    y5, y1 = [draws[i]["rect"].y0 for i in (1141, 1142)]
    colors = [(16.13, [1, 0, 0]), (16.92, [.179932, .550781, .340088])]
    rows = []
    for draw in draws:
        box = draw["rect"]
        if draw["type"] != "f" or not (249 < box.x0 < 499 and 378 < box.y0 < 475):
            continue
        for velocity, rgb in colors:
            if np.allclose(draw["fill"], rgb, atol=1e-4):
                x = (box.x0 + box.x1) / 2
                y = (box.y0 + box.y1) / 2
                mjd = (x-offset)/slope
                flux = np.exp((y1-y)/(y1-y5)*np.log(5))
                rows.append((velocity, mjd, flux))
    table = pd.DataFrame(rows, columns=["velocity_km_s", "mjd", "flux_jy"])
    if table.groupby("velocity_km_s").size().tolist() != [44, 44]:
        raise ValueError("Unexpected epoch counts")
    return table


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, help="Use a local copy of the original supplement")
    parser.add_argument("--output", type=Path, default=HERE / "torun.csv")
    args = parser.parse_args()
    if args.pdf:
        data = args.pdf.read_bytes()
    else:
        with urlopen(SOURCE, timeout=90) as response:
            data = response.read()
    extract(data).to_csv(args.output, index=False, float_format="%.5f")
