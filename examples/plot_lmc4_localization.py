"""Draw the measured LMC-4 sky offsets from the released evidence."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main(evidence, catalog, output):
    data = json.loads(evidence.read_bytes())
    ground = json.loads(catalog.read_bytes())
    fig, ax = plt.subplots(figsize=(6.4, 5), constrained_layout=True)
    rows = [r for r in data["cutouts"] if r["size"] == 13 and r["status"] == "measured"]
    normal = [r for r in rows if not r["position"]["gaussian_bound_active"]]
    failed = [r for r in rows if r["position"]["gaussian_bound_active"]]
    for subset, marker, label in [
        (normal, "o", "Individual TESS sectors"),
        (failed, "x", "Fit reached a bound"),
    ]:
        ax.scatter(
            [r["position"]["east_arcsec"] for r in subset],
            [r["position"]["north_arcsec"] for r in subset],
            marker=marker,
            s=32,
            color="#267087",
            label=label,
            zorder=3,
        )
    center = np.array(data["center_icrs_deg"])

    def tangent(ra, dec):
        return (np.array([ra, dec]) - center) * [
            3600 * np.cos(np.deg2rad(center[1])),
            3600,
        ]

    for source in data["gaia_cone"]["sources_within_5_arcsec"]:
        xy = tangent(source["ra_deg"], source["dec_deg"])
        ax.scatter(
            *xy, marker="o", facecolors="none", edgecolors="#565e66", s=35, zorder=2
        )
    ax.scatter(
        [],
        [],
        marker="o",
        facecolors="none",
        edgecolors="#565e66",
        s=35,
        label="Gaia sources within 5″ of median",
    )
    names = {
        "UCAC2   1797921": ("UCAC2 1797921", (-15, -26)),
        "OGLE LMC-LPV-20399": ("OGLE LPV", (6, 4)),
        "OGLE J050812.41-691410.1": ("Other Be candidate", (6, 4)),
    }
    for obj in ground["simbad"]["objects"]:
        if obj["main_id"] not in names:
            continue
        xy = tangent(float(obj["ra"]), float(obj["dec"]))
        label, offset = names[obj["main_id"]]
        ax.scatter(
            *xy,
            marker="*",
            s=130,
            facecolors="#e6a633",
            edgecolors="#674a18",
            linewidths=0.8,
            zorder=4,
        )
        ax.annotate(label, xy, xytext=offset, textcoords="offset points", fontsize=9)
    median = data["thirteen_pixel_summary"]["median_east_north_arcsec"]
    ax.scatter(
        *median,
        marker="+",
        color="#bd423c",
        s=110,
        linewidths=2.1,
        label="12-sector median",
        zorder=5,
    )
    ax.set(
        xlim=(-14, 16),
        ylim=(-18, 9),
        aspect="equal",
        xlabel="Offset east of original TESS centroid (arcsec)",
        ylabel="Offset north (arcsec)",
    )
    ax.grid(alpha=0.18)
    ax.legend(frameon=False, loc="upper right", fontsize=8)
    fig.savefig(output, metadata={"Date": None})
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence", type=Path)
    p.add_argument("catalog", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    main(a.evidence, a.catalog, a.output)
