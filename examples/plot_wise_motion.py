"""Plot the two blind WISE motion selections and the known-dwarf calibration."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from open_space.wise import DAYS_YEAR, REF_MJD, tangent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    fields = json.loads(args.result.read_bytes())["fields"]
    selected = [
        (field["field"], source)
        for field in fields
        for source in field["sources"]
        if source.get("disposition") in {"catalog_motion_recovered", "known_positive_control"}
    ]
    if not selected:
        raise ValueError("No confirmed control or catalog-matched track in this result")
    plt.rcParams.update(
        {
            "svg.fonttype": "none",
            "svg.hashsalt": "wise-dr3-motion-1",
            "font.size": 10,
            "font.family": "DejaVu Sans",
        }
    )
    fig, axes = plt.subplots(len(selected), 2, figsize=(9.1, 1.65 * len(selected) + 1.7), squeeze=False)
    fig.subplots_adjust(left=0.095, right=0.975, bottom=0.115, top=0.84, hspace=0.80, wspace=0.17)
    for axrow, (field, source) in zip(axes, selected):
        gaia = source["gaia"]
        center = {"ra_deg": source["ra_deg"], "dec_deg": source["dec_deg"]}
        epochs = source["epoch_measurements"]
        years = [2016 + (e["mjd"] - REF_MJD) / DAYS_YEAR for e in epochs]
        offsets = [tangent(e["ra_deg"], e["dec_deg"], center) for e in epochs]
        for index, (axis, key, label) in enumerate(
            zip(axrow, ("pmra_mas_yr", "pmdec_mas_yr"), ("RA × cos δ", "Dec"))
        ):
            axis.plot(
                years,
                [v[index] for v in offsets],
                "o",
                ms=3.4,
                color="#16536c",
                label="NEOWISE visit centroid",
            )
            interval = [min(years), max(years)]
            axis.plot(
                interval,
                [(yr - 2016) * source[key] / 1000 for yr in interval],
                color="#1f627a",
                label="NEOWISE slope",
            )
            axis.plot(
                interval,
                [(yr - 2016) * gaia[key] / 1000 for yr in interval],
                "--",
                color="#bc6339",
                label="Gaia DR3 slope",
            )
            axis.set_ylabel(label + " (arcsec)", fontsize=9)
            axis.spines[["top", "right"]].set_visible(False)
            axis.set_xticks([2014, 2016, 2018, 2020, 2022, 2024])
            axis.grid(alpha=0.18)
        title = "Known L dwarf" if field["role"] == "positive_control" else "Blind field " + field["id"]
        axrow[0].set_title(title + "    ·    Gaia " + gaia["source_id"], loc="left", fontsize=10, pad=5)
    fig.text(
        0.05, 0.958, "Infrared proper motion over a decade", fontsize=18, color="#193840", weight="bold"
    )
    fig.text(
        0.05,
        0.906,
        "Measured visit centroids and catalog slopes relative to their own positions at J2016.0.",
        fontsize=10,
        color="#41545f",
    )
    axes[-1, 0].set_xlabel("Julian year")
    axes[-1, 1].set_xlabel("Julian year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 0.002), frameon=False)
    fig.savefig(args.output, dpi=150, metadata={"Date": None} if args.output.suffix == ".svg" else None)


if __name__ == "__main__":
    main()
