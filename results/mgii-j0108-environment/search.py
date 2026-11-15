"""Search the DESI DR1 main-survey Mg II VAC. Requires numpy, scipy and astropy."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

C = 299792.458
VAC_HASH = "38c2e6281ce983a2b34b5810fdf0e3021cebaad3848cfc1cca974f185901bddf"
IDS = (39627754810908189, 39627754815095242, 39627760854893276)


def directions(ra, dec):
    a, d = np.deg2rad([ra, dec])
    return np.stack([np.cos(a) * np.cos(d), np.sin(a) * np.cos(d), np.sin(d)], axis=-1)


def selected(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != VAC_HASH:
        raise ValueError("DESI input differs from the published DR1 v1.0 checksum")
    with fits.open(path) as hdus:
        data = hdus[1].data
        # A native NumPy array keeps the input row ordering and numeric values.
        select = ((data["SURVEY"] == "main") & (data["ZWARN"] == 0)
                  & (data["LINE_SNR_MIN"] >= 3) & (data["LINE_SNR_MAX"] >= 5)
                  & (data["EW_2796"] > 0.3)
                  & (C * np.log((1 + data["Z_QSO"]) / (1 + data["Z_MGII"])) > 5000))
        return np.array(data[select])


def sky_triangles(data):
    direction = directions(data["RA"], data["DEC"])
    radius = 2 * np.sin(np.deg2rad(15 / 60) / 2)
    pairs = cKDTree(direction).query_pairs(radius, output_type="ndarray")
    distance = np.linalg.norm(direction[pairs[:, 0]] - direction[pairs[:, 1]], axis=1)
    pairs = pairs[distance > 2 * np.sin(np.deg2rad(4 / 3600) / 2)]
    graph = csr_matrix((np.ones(len(pairs), dtype=bool), pairs.T), shape=(len(data), len(data)))
    graph.sort_indices()
    triangles = []
    for i in np.flatnonzero(np.diff(graph.indptr) >= 2):
        adjacent = graph.indices[graph.indptr[i]:graph.indptr[i + 1]]
        for offset, j in enumerate(adjacent[:-1]):
            other = graph.indices[graph.indptr[j]:graph.indptr[j + 1]]
            common = np.intersect1d(adjacent[offset + 1:], other, assume_unique=True)
            triangles.extend((i, j, k) for k in common)
    triangles = np.asarray(triangles)
    length = [np.linalg.norm(direction[triangles[:, i]] - direction[triangles[:, (i+1) % 3]],
                             axis=1) for i in range(3)]
    diam = np.rad2deg(2 * np.arcsin(np.max(length, axis=0) / 2)) * 60
    return triangles, diam


def exchange_blocks(data, regional):
    # Keep all redshifts and EWs on a given donor spectrum together and preserve multiplicity.
    _, obj = np.unique(data["TARGETID"], return_inverse=True)
    ordered = np.lexsort((data["Z_MGII"], obj))
    _, start, count = np.unique(obj[ordered], return_index=True, return_counts=True)
    base = ordered[start]
    zq = data["Z_QSO"][base]
    region = np.zeros(len(zq), dtype=int)
    if regional:
        region = ((data["RA"][base] % 360 // 60).astype(int)
                  + 6 * ((data["DEC"][base] + 90) // 30).astype(int))
    qbin = np.floor(zq * 10).astype(int)
    labels = np.lexsort((zq, qbin, region, count))
    boundaries = np.r_[0, np.flatnonzero((np.diff(count[labels]) != 0)
                        | (np.diff(region[labels]) != 0) | (np.diff(qbin[labels]) != 0)) + 1,
                        len(labels)]
    for lo, hi in zip(boundaries[:-1], boundaries[1:]):
        for group in np.array_split(labels[lo:hi], (hi - lo + 79) // 80):
            sources = ordered[start[group, None] + np.arange(count[group[0]])]
            vmax_z = np.max(data["Z_MGII"][sources], axis=1)
            safe = C * np.log((1 + np.min(zq[group])) / (1 + vmax_z)) > 5000
            if np.any(safe):
                yield sources[safe]


def shuffle_test(data, tri, n_sim, regional):
    z, ew = np.log1p(data["Z_MGII"]), data["EW_2796"]
    used = np.zeros(len(data), dtype=bool)
    used[tri] = True
    blocks = []
    for group in exchange_blocks(data, regional):
        needed = group[used[group].any(axis=1)]
        if len(needed):
            blocks.append((needed, group))
    rng = np.random.default_rng(8871)
    source = np.arange(len(data))
    hits = 0
    total = 0
    for _ in range(n_sim):
        for needed, donors in blocks:
            source[needed] = donors[rng.permutation(len(donors))[:len(needed)]]
        trial = source[tri]
        n = int(np.count_nonzero((C * np.ptp(z[trial], axis=1) <= 250)
                                 & (np.min(ew[trial], axis=1) >= 1)))
        hits += n > 0
        total += n
    return {"simulations": n_sim, "seed": 8871, "with_one_or_more": hits,
            "total_triplets": total, "plus_one_tail_fraction": (hits+1)/(n_sim+1)}


def write_csv(path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("catalog", type=Path)
    parser.add_argument("--out", type=Path, default=Path("."))
    parser.add_argument("--shuffles", type=int, default=5000)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    data = selected(args.catalog)
    tri, diam = sky_triangles(data)
    ew = np.min(data["EW_2796"][tri], axis=1)
    vel = C * np.ptp(np.log1p(data["Z_MGII"][tri]), axis=1)
    inds = np.flatnonzero((ew >= 1) & (vel <= 500))
    ordered = inds[np.argsort(diam[inds])]
    write_csv(args.out / "all-strong-triplets.csv",
              ["diameter_arcmin", "velocity_span_kms", "minimum_EW2796_A",
               "targetid_1", "targetid_2", "targetid_3"],
              ((diam[j], vel[j], ew[j], *data["TARGETID"][tri[j]]) for j in ordered))
    cols = ["TARGETID", "RA", "DEC", "Z_QSO", "Z_MGII", "Z_MGII_ERR_LOW", "Z_MGII_ERR_HIGH",
            "EW_2796", "EW_2796_ERR_LOW", "EW_2796_ERR_HIGH", "EW_2803",
            "EW_2803_ERR_LOW", "EW_2803_ERR_HIGH", "LINE_SNR_MIN", "LINE_SNR_MAX"]
    j0108 = [np.flatnonzero((data["TARGETID"] == i) & (abs(data["Z_MGII"] - .96) < .002)).item()
             for i in IDS]
    write_csv(args.out / "j0108-absorbers.csv", cols, ([data[k][j] for k in cols] for j in j0108))
    summary = {"input_sha256": VAC_HASH, "systems": len(data),
               "sightlines": len(np.unique(data["TARGETID"])), "triangles_15arcmin": len(tri),
               "triangles_6arcmin": int(sum(diam <= 6)), "strong_15arcmin_500kms": len(inds),
               "strong_6arcmin_250kms": int(sum((diam <= 6) & (ew >= 1) & (vel <= 250)))}
    if args.shuffles:
        for label, regional in [("all_sky_shuffle", False), ("regional_shuffle", True)]:
            summary[label] = shuffle_test(data, tri[diam <= 6], args.shuffles, regional)
    (args.out / "selection.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
