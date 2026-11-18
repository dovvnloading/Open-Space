"""Exploratory scale sensitivity on the frozen three-cell footprints."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

from screen import HERE, load_current, locate, values


WIDTHS_KM = (25, 40, 60)


def compare_scales(root, frozen, sources):
    grid = {}
    rows = []
    by_date = {entry["date"]: entry for entry in sources if entry["sensor"] == "amsr3"}
    for day in frozen["daily"]:
        openings = day["variants"]["core3"]["openings"]
        if not openings:
            continue
        entry = by_date[day["date"]]
        frame = load_current(locate(root, entry, False), entry, grid)
        eligible = frame.valid & frame.offshore
        quantity = np.where(eligible, frame.effective, 0)
        smoothed = {}
        for width in WIDTHS_KM:
            sigma = width / (10 * np.sqrt(8 * np.log(2)))
            total = gaussian_filter(quantity, sigma)
            weight = gaussian_filter(eligible.astype(float), sigma)
            smoothed[width] = np.divide(total, weight, out=np.full_like(total, np.nan), where=weight >= 0.8)
        for opening in openings:
            meds = {}
            for width, image in smoothed.items():
                core = values(image, opening["core"])
                meds[str(width)] = round(float(np.median(core)), 3) if np.isfinite(core).all() else None
            rows.append({
                "id": opening["id"],
                "center_deg": opening["center_deg"],
                "amsr3_unsmoothed_median_pct": opening["core_median_pct"],
                "amsr3_smoothed_medians_pct": meds,
                "amsr2_median_pct": opening["cross_sensor"]["amsr2"].get("core_median_pct"),
                "ssmis_median_pct": opening["cross_sensor"]["ssmis"].get("core_median_pct"),
            })
    return rows


def plot_chain(rows, filename):
    import matplotlib.pyplot as plt

    days = ["15 Nov", "16 Nov", "17 Nov"]
    ids = ("2026-11-15/1", "2026-11-16/1", "2026-11-17/1")
    by_id = {item["id"]: item for item in rows}
    selected = [by_id[item] for item in ids]
    fig, ax = plt.subplots(figsize=(7.7, 4.1), layout="constrained")
    ax.plot(days, [d["amsr3_unsmoothed_median_pct"] for d in selected], "o-",
            color="#126c6d", label="AMSR3")
    ax.plot(days, [d["amsr2_median_pct"] for d in selected], "s--",
            color="#2a4f85", label="AMSR2 (17 Nov unavailable)")
    ax.plot(days, [d["ssmis_median_pct"] for d in selected], "D-",
            color="#9c2a45", label="SSMIS")
    for width, alpha in (("40", 0.5), ("60", 0.8)):
        ax.plot(days, [d["amsr3_smoothed_medians_pct"][width] for d in selected], ":",
                color="#73552d", alpha=alpha, linewidth=2.5, label=f"AMSR3 with {width} km smoothing")
    ax.set_title("The same Antarctic sea-ice patch looks different across sensors", loc="left", fontsize=13)
    ax.set_ylabel("Median estimated ice concentration (%)")
    ax.set_ylim(0, 100)
    ax.grid(axis="y", color="#e3e8eb")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower left", ncols=2, fontsize=8, frameon=False)
    fig.savefig(filename, dpi=140)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--screen", type=Path, default=HERE / "screen-results.json")
    args = parser.parse_args()
    frozen = json.loads(args.screen.read_text(encoding="utf-8"))
    sources = json.loads((HERE / "sources.json").read_text(encoding="utf-8"))
    rows = compare_scales(args.data, frozen, sources)
    args.output.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8", newline="\n")
    if args.plot:
        plot_chain(rows, args.plot)


if __name__ == "__main__":
    main()
