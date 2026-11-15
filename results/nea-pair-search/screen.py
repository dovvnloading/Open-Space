"""Rank pairs in a JPL SBDB near-Earth asteroid snapshot. Requires NumPy and SciPy."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.spatial import cKDTree

FIELDS = "spkid,pdes,full_name,epoch,e,a,q,i,om,w,ma,n,condition_code,data_arc,n_obs_used,first_obs,last_obs,H,orbit_id"
QUERY = "https://ssd-api.jpl.nasa.gov/sbdb_query.api?" + urlencode({
    "fields": FIELDS, "sb-group": "neo", "sb-kind": "a", "full-prec": "true",
})
THRESHOLD = 0.025


def orbit_vectors(table):
    def column(key):
        return np.array([float(row[key]) for row in table])

    inc, node, peri = np.radians([column("i"), column("om"), column("w")])
    ci, si = np.cos(inc), np.sin(inc)
    cn, sn = np.cos(node), np.sin(node)
    cp, sp = np.cos(peri), np.sin(peri)
    pole = np.column_stack((si * sn, -si * cn, ci))
    pericenter = np.column_stack((cn * cp - sn * sp * ci,
                                  sn * cp + cn * sp * ci, sp * si))
    return np.column_stack((column("q"), pole, column("e")[:, None] * pericenter))


def rank(raw, output):
    catalog = json.loads(raw)
    rows = [dict(zip(catalog["fields"], values)) for values in catalog["data"]]
    vectors = orbit_vectors(rows)
    indices = cKDTree(vectors).query_pairs(THRESHOLD, output_type="ndarray")
    scores = np.linalg.norm(vectors[indices[:, 0]] - vectors[indices[:, 1]], axis=1)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["rank", "score", "first", "second", "epoch1_jd_tdb",
                         "epoch2_jd_tdb", "arc1_days", "arc2_days", "U1", "U2", "H1", "H2"])
        for rank_number, k in enumerate(np.argsort(scores), start=1):
            first, second = (rows[index] for index in indices[k])
            writer.writerow([rank_number, f"{scores[k]:.9f}", first["pdes"], second["pdes"],
                             first["epoch"], second["epoch"], first["data_arc"], second["data_arc"],
                             first["condition_code"], second["condition_code"], first["H"], second["H"]])
    n = len(rows)
    print(json.dumps({"objects": n, "possible_pairs": n * (n - 1) // 2,
                      "pairs_below_threshold": len(scores), "input_sha256": hashlib.sha256(raw).hexdigest(),
                      "query": QUERY}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path, help="raw JSON snapshot or cache destination")
    parser.add_argument("output", type=Path, help="destination CSV")
    args = parser.parse_args()
    if args.snapshot.exists():
        data = args.snapshot.read_bytes()
    else:
        with urlopen(QUERY, timeout=180) as response:
            data = response.read()
        json.loads(data)
        args.snapshot.write_bytes(data)
    rank(data, args.output)
