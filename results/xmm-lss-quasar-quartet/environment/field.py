"""Compare the quasar field to nearby apertures in the Chen et al. photo-z catalog."""
import argparse
import csv
import io
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.spatial import cKDTree

RA, DEC, Z = 36.31513, -4.78245, 1.9295
MEMBERS = Path(__file__).resolve().parents[1] / "members.csv"
URL = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?" + urlencode({
    "-source": "J/MNRAS/478/2132/tablec", "-c": f"{RA} {DEC}", "-c.rm": "40",
    "-out": "VIDEOID,RAJ2000,DEJ2000,zph,q_zph,zsp", "-out.max": "unlimited",
})


def coordinates(ra, dec):
    return np.column_stack(((ra - RA) * 60 * np.cos(np.deg2rad(DEC)), (dec - DEC) * 60))


def read_catalog(path):
    if not path.exists():
        with urlopen(URL, timeout=180) as response:
            content = response.read()
        path.write_bytes(content)
    lines = [line for line in path.read_text().splitlines() if line and not line.startswith("#")]
    rows = csv.DictReader(io.StringIO("\n".join([lines[0], *lines[3:]])), delimiter="\t")
    return np.array([(float(r["RAJ2000"]), float(r["DEJ2000"]), float(r["zph"]))
                     for r in rows if r["zph"].strip()])


def compare(data, member_file):
    ra, dec, photo_z = data.T
    points = coordinates(ra, dec)
    tree = cKDTree(points)
    with member_file.open() as handle:
        members = list(csv.DictReader(handle))
    member_xy = coordinates(np.array([float(row["ra_deg"]) for row in members]),
                            np.array([float(row["dec_deg"]) for row in members]))
    # Remove catalog detections within 2 arcsec of a quartet member.
    distinct = cKDTree(member_xy).query(points)[0] > 2 / 60
    control = np.stack(np.meshgrid(np.arange(-34, 35, 2), np.arange(-34, 35, 2)),
                       axis=-1).reshape(-1, 2)
    distance = np.linalg.norm(control, axis=1)
    control = control[(distance < 34) & (distance > 5)]
    rows = []
    for radius in (1, 2, 3):
        target = np.array(tree.query_ball_point([0, 0], radius), dtype=int)
        target = target[distinct[target]]
        pools = tree.query_ball_point(control, radius)
        counts = np.array([len(p) for p in pools])
        matched = (counts >= .8 * target.size) & (counts <= 1.2 * target.size) & (counts > 0)
        for width in (.05, .10, .15):
            selected = np.abs(photo_z - Z) <= width
            count = int(selected[target].sum())
            fractions = np.array([selected[p].sum() for p in pools])[matched] / counts[matched]
            rows.append((radius, width, target.size, count, count / target.size,
                         len(fractions), np.median(fractions), np.mean(fractions >= count / target.size)))
    return rows, int((~distinct).sum())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=Path("chen-table-c.tsv"))
    parser.add_argument("--members", type=Path, default=MEMBERS)
    parser.add_argument("--output", type=Path, default=Path("apertures.csv"))
    args = parser.parse_args()
    data = read_catalog(args.table)
    rows, removed = compare(data, args.members)
    with args.output.open("w", newline="") as handle:
        out = csv.writer(handle)
        out.writerow(["radius_arcmin", "dz_half_width", "field_n", "field_in_slice",
                      "field_fraction", "control_apertures", "median_control_fraction",
                      "fraction_of_controls_at_least_as_high"])
        out.writerows(rows)
    print(f"{len(data)} input rows, {removed} counterparts excluded")
