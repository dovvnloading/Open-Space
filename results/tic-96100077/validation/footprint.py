"""Map each 9x9 TESSCut to physical detector columns using its FITS headers."""

import argparse
from pathlib import Path

import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parent


def physical(header, axis, pixel):
    return (
        header[f"{axis}CRV5P"]
        + (pixel - header[f"{axis}CRP5P"]) * header[f"{axis}CDL5P"]
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=ROOT.parent / "raw")
    args = parser.parse_args()
    orbit = pd.read_csv(ROOT.parent / "data/gaia-orbit.csv").iloc[0]
    results = []
    for sector in pd.read_csv(ROOT.parent / "outputs/sectors.csv").sector:
        with fits.open(args.raw_dir / f"s{sector:03}.fits") as hdus:
            header = hdus[1].header
            sky = hdus[2].header
            x, y = WCS(sky).all_world2pix(orbit.ra, orbit.dec, 1)
            col, row = physical(header, 1, x), physical(header, 2, y)
            results.append(
                dict(
                    sector=sector,
                    first_col=physical(header, 1, 1),
                    last_col=physical(header, 1, sky["NAXIS1"]),
                    target_col=col,
                    target_row=row,
                    in_science_area=45 <= col <= 2092 and 1 <= row <= 2048,
                )
            )
    pd.DataFrame(results).to_csv(
        ROOT / "footprint.csv", index=False, float_format="%.4f"
    )
