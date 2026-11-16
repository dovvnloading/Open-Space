"""Screen the detected DR1 population. Write every eligible row, including ordinary spectra."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.spatial import cKDTree

from comparison.match_catalogues import unit_xyz


def run(release, matches, output):
    cat = Table.read(release / "sources.fits").to_pandas()
    cat["S_Code"] = cat.S_Code.str.decode("ascii").str.strip()
    fields = pd.read_csv(release / "source_fields.csv", index_col="row")
    cat = cat.join(fields)
    for prefix in ("lotss", "nvss", "wenss", "gmrtas150m", "vlssr", "rrs8c38mhz", "gleamegcat"):
        other = pd.read_csv(matches / f"match_{prefix}.csv.gz", index_col="row")
        if not other.index.equals(cat.index):
            raise ValueError("External catalogue row order differs from the native FITS")
        cat = cat.join(other)
    positions = unit_xyz(cat.RA, cat.DEC)
    chord = cKDTree(positions).query(positions, k=2)[0][:, 1]
    cat["neighbour_arcsec"] = np.degrees(2 * np.arcsin(chord / 2)) * 3600
    cuts = [
        ("all_rows", np.full(len(cat), True)), ("good_field", cat.quality == "good"),
        ("within_3p5_deg", cat.separation_deg < 3.5),
        ("simple_bright_isolated", (cat.S_Code == "S") & (cat.Peak_flux / cat.Isl_rms > 7)
         & (cat.Total_flux / cat.Peak_flux < 1.7) & (cat.neighbour_arcsec > 90)),
        ("compact_144", (cat.lotss_sep_as < 15) & (cat.lotss_flux_jy / cat.lotss_sum75_jy > .9)
         & (cat.lotss_dc_maj < 12) & (cat.lotss_flux_jy / cat.lotss_peak_jy_beam < 1.45)),
        ("unique_1400", (cat.nvss_n_match == 1) & (cat.nvss_sep_as < 18)),
    ]
    eligible = np.full(len(cat), True)
    counts = []
    for cut, mask in cuts:
        eligible &= mask
        counts.append((cut, int(eligible.sum())))
    frame = cat[eligible].copy()
    alpha = np.log(frame.nvss_flux_jy / frame.lotss_flux_jy) / np.log(1400 / 144)
    # A screening contrast; formal Gaussian sigmas cannot capture the field response.
    frame["alpha_144_1400"] = alpha
    frame["pred_over_23"] = frame.lotss_flux_jy * (23 / 144)**alpha / frame.Total_flux
    frame["pred_over_74"] = frame.lotss_flux_jy * (73.8 / 144)**alpha / frame.vlssr_flux_jy
    keep = ["RA", "DEC", "field", "separation_deg", "Total_flux", "E_Total_flux", "Peak_flux",
            "Isl_rms", "neighbour_arcsec", "lotss_name", "lotss_dc_maj", "lotss_sum75_jy",
            "alpha_144_1400", "pred_over_23", "pred_over_74"]
    keep += [col for col in frame if col.endswith(("_flux_jy", "_err_jy", "_sep_as"))]
    output.mkdir(parents=True, exist_ok=True)
    frame[keep].to_csv(output / "detected_sample.csv", index_label="row", float_format="%.8g")
    pd.DataFrame(counts, columns=["cut", "remaining"]).to_csv(output / "selection.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ["release_dir", "match_dir", "output"]:
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    run(args.release_dir, args.match_dir, args.output)
