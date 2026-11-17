"""Propagate the published-centroid grid to each unique NOIRLab file-search midpoint.

Duplicate reductions of the same exposure are grouped. A Data Lab four-corner polygon
match, when available, is geometrical coverage and does not account for bad pixels.
"""

from pathlib import Path
import csv, json
from collections import defaultdict
from matplotlib.path import Path as Polygon
import numpy as np
from decam_difference import p, HERE


def main(src, dest):
    states = p.family(
        json.loads(Path(p.__file__).with_name("orbit_input.json").read_text())
    )
    files = json.loads((src / "all_decam_files.response").read_text())[1:]
    idx = defaultdict(list)
    obs = defaultdict(list)
    for c in csv.DictReader((src / "nsa_region_valid.response").open()):
        idx[c["fileref"]].append(c)
    for f in files:
        obs[f["dateobs_center"]].append(f)
    with (dest / "decam_epochs.csv").open("w", newline="") as fp:
        wr = csv.writer(fp)
        wr.writerow(
            [
                "midpoint_utc",
                "filter",
                "exposure_s",
                "proposal",
                "object",
                "release_date",
                "archive_products",
                "ra_min",
                "ra_max",
                "dec_min",
                "dec_max",
                "boresight_ra",
                "boresight_dec",
                "n_grid_nodes",
                "legacy_footprint_rows",
                "nodes_in_legacy_CCD_polygons",
                "min_boresight_distance_deg",
            ]
        )
        for date, rows in sorted(obs.items()):
            row = rows[0]
            a, d = p.sky(states, p.epoch(date))
            pts = np.stack((a, d), axis=1)
            ccds = [c for r in rows for c in idx[Path(r["archive_filename"]).name]]
            hits = np.zeros(len(a), dtype=bool)
            for c in ccds:
                corners = [
                    [float(c[f"ra{i}"]), float(c[f"dec{i}"])] for i in (1, 2, 3, 4)
                ]
                hits |= Polygon(corners).contains_points(pts)
            center = p.unit(row["ra_center"], row["dec_center"])
            s = p.direction(states, p.epoch(date))
            angles = np.rad2deg(2 * np.arcsin(p.norm(s - center) / 2))
            wr.writerow(
                [
                    date,
                    row["ifilter"].split()[0],
                    row["exposure"],
                    row["proposal"],
                    row["OBJECT"],
                    row["release_date"],
                    len(rows),
                    *[round(float(v), 6) for v in [a.min(), a.max(), d.min(), d.max()]],
                    row["ra_center"],
                    row["dec_center"],
                    len(a),
                    len(ccds),
                    int(hits.sum()) if ccds else "",
                    round(angles.min(), 5),
                ]
            )


if __name__ == "__main__":
    main(HERE / "receipts", HERE / "reduced")
