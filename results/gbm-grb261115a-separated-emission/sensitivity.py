"""Count-space sensitivity of the frozen CTIME coincidence test."""

from argparse import ArgumentParser
import csv
import json
from pathlib import Path

import numpy as np
from scipy.special import gammaincinv

from screen import exact_tail


def minimum_count(exposure, sides, alpha):
    """Smallest integer count for which both one-sided reference tests pass."""
    low, high = 0, 1

    def exceeds_threshold(count):
        return max(exact_tail(count, exposure, side) for side in sides.values()) > alpha

    while exceeds_threshold(high):
        high *= 2
    while high - low > 1:
        middle = (low + high) // 2
        if exceeds_threshold(middle):
            low = middle
        else:
            high = middle
    return high


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("screen_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    source = json.loads((args.screen_dir / "screen.json").read_text(encoding="utf-8"))
    levels = {"n9": source["unadjusted_joint_threshold"],
              "n6": source["unadjusted_joint_threshold"], "n7": source["family_alpha"]}
    # A union bound applies to the three failures regardless of detector dependence.
    required_marginal_power = 1 - 0.10 / len(levels)
    groups = {}
    with (args.screen_dir / "windows.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["region"] not in ("pre", "late"):
                continue
            additions = []
            for detector, alpha in levels.items():
                exposure = float(row[f"{detector}_exposure_s"])
                sides = source["backgrounds"][detector]
                count = minimum_count(exposure, sides, alpha)
                # The lower measured sideband is the fixed count-model baseline.
                background_mean = exposure * min(side["rate"] for side in sides.values())
                required_mean = gammaincinv(count, required_marginal_power)
                additions.append(max(0.0, float(required_mean) - background_mean))
            key = (row["region"], float(row["nominal_width_s"]))
            groups.setdefault(key, []).append((max(additions), row["direction_passes"] == "True"))

    result = {
        "model": "Fixed observed sidebands; Poisson on-window counts at lower sideband rate plus equal added mean counts in each of n9, n6 and n7. Directional veto omitted from modeled power.",
        "min_joint_power": 0.90,
        "min_marginal_power": required_marginal_power,
        "groups": [{
            "region": region, "nominal_width_s": width, "windows": len(rows),
            "observed_direction_clear": sum(clear for _, clear in rows),
            "added_mean_counts_per_detector_min_median_max": [
                float(value) for value in np.quantile([added for added, _ in rows], [0, 0.5, 1])
            ],
        } for (region, width), rows in sorted(groups.items())],
    }
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
