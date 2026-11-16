"""Count visits and test apparent white-dwarf flux deficits against field controls."""
import csv
import gzip
from collections import defaultdict
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent


def read(name):
    path = ROOT/name
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream))


matches = read("matches.csv")
detections = read("detection.csv.gz")
forced = read("forced_photometry.csv.gz")
print(f"{len(matches)} DIA objects, {len({r['gaia'] for r in matches})} white dwarfs, "
      f"{len(detections)} detections, {len(forced)} forced fluxes")
groups = defaultdict(list)
for row in read("control_visits.csv"):
    groups[(row["target"], row["visit"], row["detector"])].append(row)
for (target, visit, detector), rows in sorted(groups.items()):
    t = next(r for r in rows if r["object"] == target)
    comp = [float(r["depth"]) for r in rows if r["object"] != target]
    print(visit, detector, f"target depth {float(t['depth']):.3f}",
          f"neighbors >50% low {sum(x>.5 for x in comp)}/{len(comp)}",
          f"neighbor median deficit {median(comp):.3f}",
          f"projected CCD pixel ({t['x_est']}, {t['y_est']})")
