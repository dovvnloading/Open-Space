"""Rebuild the frozen orbit-neighbor table from the named MPCORB snapshot."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

DATA = Path(__file__).resolve().parent / "data"


def read_orbits(path):
    names, elements = [], []
    with gzip.open(path, "rt", encoding="latin1") as stream:
        for line in stream:
            if len(line) < 200 or line[20:25] != "K2669":
                continue
            a = float(line[92:103])
            if int(line[123:126]) < 2 or not 2 < a < 5.6:
                continue
            names.append(line[166:194].strip())
            elements.append(
                (
                    a,
                    float(line[70:79]),
                    float(line[59:68]),
                    float(line[48:57]),
                    float(line[37:46]),
                    float(line[26:35]),
                )
            )
    return names, np.array(elements).T


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    metadata = json.loads((DATA / "sources.json").read_text())["MPCORB"]
    if (
        hashlib.sha256(args.snapshot.read_bytes()).hexdigest()
        != metadata["gzip_sha256"]
    ):
        parser.error(
            "Historical replay requires the MPCORB gzip identified in data/sources.json."
        )
    names, (a, e, inc, node, peri, anomaly) = read_orbits(args.snapshot)
    i, om, w = np.deg2rad([inc, node, peri])
    xy = np.column_stack(
        (
            np.sqrt(1.25) * np.log(a),
            np.sqrt(2) * e * np.cos(om + w),
            np.sqrt(2) * e * np.sin(om + w),
            2 * np.sqrt(2) * np.sin(i / 2) * np.cos(om),
            2 * np.sqrt(2) * np.sin(i / 2) * np.sin(om),
        )
    )
    pairs = cKDTree(xy).query_pairs(0.0014, output_type="ndarray")
    left, right = pairs.T
    metric = (
        np.linalg.norm(xy[left] - xy[right], axis=1)
        * 29785
        / np.sqrt((a[left] + a[right]) / 2)
    )
    longitude = node + peri + anomaly
    delta = abs((longitude[left] - longitude[right] + 180) % 360 - 180)
    with args.output.open("w") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            ["first", "second", "metric_m_s", "delta_longitude_deg", "delta_a_au"]
        )
        for k in np.argsort(metric):
            l, r = pairs[k]
            writer.writerow(
                [
                    names[l],
                    names[r],
                    format(metric[k], ".10g"),
                    format(delta[k], ".10g"),
                    format(a[l] - a[r], ".10g"),
                ]
            )
    print(f"{len(names)} eligible orbits, {len(pairs)} retained pairs.")


if __name__ == "__main__":
    main()
