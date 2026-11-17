"""Plot the actual 60-micron matched-filter amplitudes against detector offset."""
import csv
import json
import re
import tarfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from scan_tracks import filtered

ROOT = Path(__file__).resolve().parent


def main():
    root = ROOT / "scanpi"
    with (root / "track_measurements.csv").open() as stream:
        values = {int(r["scan"]): r for r in csv.DictReader(stream) if r["band_um"] == "60"}
    table = Table.read(root / "b3_summary.tbl", format="ascii.ipac")
    grid = np.linspace(-8, 5, 781)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.9), sharey=True)
    blank = json.loads((root / "track_summary.json").read_text())
    with tarfile.open(root / "profiles.tar.gz") as data:
        for ax, group, title in zip(axes, ("July", "Aug_Sep"), ("July 9–10", "August 28 and September 1")):
            profiles, weights, ranges = [], [], []
            for exposure in table:
                record = values.get(int(exposure["scan"]))
                if record is None:
                    continue
                month = int(record["utc"][5:7])
                if (group == "July" and month != 7) or (group == "Aug_Sep" and month < 8):
                    continue
                filename = Path(re.search(r'href="([^"]+)"', exposure["table"])[1]).name
                lines = data.extractfile("b3/" + filename).read().decode().splitlines()
                x, y = np.array([line.split()[:2] for line in lines if line.startswith(" ")], dtype=float).T
                profiles.append(np.interp(grid, x, filtered(x, y, 1.5)))
                weights.append(1 / max(0.04, float(record["control_robust_scatter_Jy"]))**2)
                ranges.extend(float(record[key]) for key in ("july_gated_along_min_arcmin", "july_gated_along_max_arcmin"))
            curve = np.average(profiles, axis=0, weights=weights)
            central = np.interp(0., grid, curve)
            assert np.isclose(central, blank["3_"+group+"_sign+1"]["fixed_Jy"])
            ax.axhline(0, lw=.7, color="#637582")
            ax.axvspan(min(ranges), max(ranges), facecolor="#f1b46c", alpha=.45, label="Predicted moving positions")
            ax.axvline(0, color="#476475", lw=.9, ls="--", label="Catalogued IRAS position")
            for y in profiles:
                ax.plot(grid, y, lw=.7, color="#667d87", alpha=.33)
            ax.plot(grid, curve, color="#005e71", lw=2.1, label="Weighted scan profile")
            ax.scatter([0], [central], s=25, color="#005e71", zorder=4)
            ax.annotate(f"{central:.3f} Jy", (0, central), xytext=(5, 9), textcoords="offset points", fontsize=9, color="#005e71")
            ax.set_title(f"{title}  ·  {len(profiles)} crossings", fontsize=11, loc="left")
            ax.set_xlim(grid[0], grid[-1])
            ax.spines[["top", "right"]].set_visible(False)
            ax.set_xlabel("Along-scan offset (arcmin)")
    axes[0].set_ylabel("Fitted point-source flux (Jy)")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.suptitle("IRAS measured the late 60 µm signal at the original sky position", x=.075, ha="left", fontsize=14)
    fig.subplots_adjust(left=.075, right=.99, bottom=.26, top=.79, wspace=.13)
    fig.savefig(ROOT/"scan_profiles.png", dpi=155)
    plt.close(fig)


if __name__ == "__main__":
    main()
