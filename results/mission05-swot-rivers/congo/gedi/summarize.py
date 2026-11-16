"""Recompute quality exclusion and footprint-level diagnostics from the archived extract."""
import csv
import gzip
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent


def main():
    by_group = defaultdict(list)
    accepted = []
    with gzip.open(ROOT / "gedi-corridor-shots.csv.gz", "rt", newline="") as f:
        for row in csv.DictReader(f):
            if float(row["offset_m"]) > 450:
                continue
            by_group[row["utc"][:10], row["reach_id"]].append(row)
    data = []
    for (date, reach), rows in sorted(by_group.items()):
        water = [r for r in rows if r["jan_scl"] == r["jul_scl"] == "6"]
        strict = [r for r in water if int(r["flags"]) == 130]
        interior = [r for r in strict if r["jan_water_3x3_count"] == r["jul_water_3x3_count"] == "9"]
        accepted.extend(interior)
        values = [float(r["elev_lowestmode_m"]) for r in interior]
        data.append([date, reach, len(rows), len(water), len(strict), len(interior),
                     len({r["beam"] for r in interior}),
                     round(min(values), 3) if values else "", round(median(values), 3) if values else "",
                     round(max(values), 3) if values else ""])
    with (ROOT / "screen.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date_utc", "reach_id", "n_within_450m", "n_both_dates_water",
                    "n_both_water_quality_ok", "n_interior_water_quality_ok", "n_distinct_beams",
                    "min_ellipsoidal_m", "median_ellipsoidal_m", "max_ellipsoidal_m"])
        w.writerows(data)
    grouped = defaultdict(list)
    for r in accepted:
        grouped[r["utc"][:10], r["reach_id"], r["beam"]].append(r)
    points = []
    for (date, reach, beam), shots in sorted(grouped.items()):
        points.append([date, reach, beam, len(shots),
                       round(median(float(r["along_km_from_downstream_end"]) for r in shots), 5),
                       round(median(float(r["elev_lowestmode_m"]) for r in shots), 4)])
    with (ROOT / "accepted-beam-medians.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date_utc", "reach_id", "beam", "n_shots", "median_along_km", "median_ellipsoidal_m"])
        w.writerows(points)
    sample = [r for r in points if r[0] == "2019-07-20" and r[1] == "13217000061"]
    secants = [(b[5] - a[5]) / (b[4] - a[4]) for a,b in combinations(sample, 2)
               if not math.isclose(a[4], b[4])]
    print("2019-07-20 downstream reach", len(sample), "beams", sum(r[3] for r in sample),
          "shots. Median pairwise ellipsoid slope m/km:", median(secants) if secants else None)
    for row in data:
        if row[5]:
            print(*row)


if __name__ == "__main__":
    main()
