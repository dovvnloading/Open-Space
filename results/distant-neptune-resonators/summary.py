"""Rebuild the current orbit comparison from the bundled per-clone evidence."""

from pathlib import Path
import csv
from collections import defaultdict

ROOT = Path(__file__).resolve().parent
WORK = ROOT / "work"


def summarize():
    records = defaultdict(dict)
    with (WORK / "mpc-clones.csv").open() as file:
        for r in csv.DictReader(file):
            if r["method"] != "whfast" or int(r["draws"]) == 0:
                continue
            key = (
                r["object"],
                int(r["draws"]),
                abs(float(r["step_yr"])),
                abs(int(r["years"])),
            )
            records[key][(int(r["clone"]), int(r["years"]) > 0)] = r
    rows = []
    for (name, n, dt, t), data in sorted(records.items()):
        spans = []
        for i in range(n + 1):
            pair = [data[i, s] for s in (False, True)]
            span = max(float(x["relative_max_deg"]) for x in pair) - min(
                float(x["relative_min_deg"]) for x in pair
            )
            spans.append(span)
        rows.append(
            [
                name,
                n,
                dt,
                t,
                round(spans[0], 4),
                sum(s < 355 for s in spans[1:]),
                sum(s < 360 for s in spans[1:]),
            ]
        )
    return rows


if __name__ == "__main__":
    with (WORK / "current-summary.csv").open("w") as f:
        wr = csv.writer(f)
        wr.writerow(
            [
                "object",
                "draws",
                "step_yr",
                "half_span_yr",
                "nominal_joint_range_deg",
                "joint_below_355",
                "joint_below_360",
            ]
        )
        wr.writerows(summarize())
