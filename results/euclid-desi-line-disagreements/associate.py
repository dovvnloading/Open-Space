"""Crossmatch DESI development targets to complete downloaded Euclid MER catalogs."""

from argparse import ArgumentParser
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

from astropy.io import fits
import numpy as np
from scipy.spatial import cKDTree


def sphere(ra, dec):
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack((np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)))


def read_mer(entries, root, bounds):
    catalog = {}
    tiles = {}
    left, right, bottom, top = bounds
    for item in entries:
        path = root / item["file"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"MER digest mismatch for {path.name}")
        with fits.open(path, memmap=True) as hdus:
            table = next(hdu.data for hdu in hdus if isinstance(hdu, fits.BinTableHDU)
                         and {"OBJECT_ID", "RIGHT_ASCENSION", "DECLINATION"} <= set(hdu.columns.names))
            select = ((table["RIGHT_ASCENSION"] >= left) & (table["RIGHT_ASCENSION"] < right)
                      & (table["DECLINATION"] >= bottom) & (table["DECLINATION"] < top))
            for row in table[select]:
                ident = int(row["OBJECT_ID"])
                values = (float(row["RIGHT_ASCENSION"]), float(row["DECLINATION"]))
                if ident in catalog and catalog[ident] != values:
                    raise ValueError(f"Reused Euclid ID at different coordinates: {ident}")
                catalog[ident] = values
                tiles.setdefault(ident, set()).add(item["tile"])
    return catalog, tiles


def associate(ground, mer, tiles):
    ids = sorted(mer)
    sky = np.asarray([mer[ident] for ident in ids])
    m_vec = sphere(sky[:, 0], sky[:, 1])
    g_vec = sphere([float(row["mean_fiber_ra"]) for row in ground],
                   [float(row["mean_fiber_dec"]) for row in ground])
    radius = 2 * np.sin(np.deg2rad(3 / 3600) / 2)
    groups = cKDTree(m_vec).query_ball_point(g_vec, radius)
    edges = []
    for i, group in enumerate(groups):
        neighbors = []
        for j in group:
            chord = np.linalg.norm(g_vec[i] - m_vec[j])
            separation = 2 * np.rad2deg(np.arcsin(chord / 2)) * 3600
            neighbors.append((separation, ids[j]))
        neighbors.sort()
        edges.append(neighbors)
    ground_per_mer = Counter(ident for group in edges for sep, ident in group if sep <= 1)

    result = []
    for row, neighbors in zip(ground, edges):
        close = [(sep, ident) for sep, ident in neighbors if sep <= 1]
        if len(close) != 1:
            preferred = None
        else:
            sep, ident = close[0]
            preferred = ident if sep <= 0.5 and ground_per_mer[ident] == 1 else None
        valid_quality = (
            float(row["deltachi2"]) >= 40 and int(row["coadd_fiberstatus"]) == 0
            and np.isfinite(float(row["zerr"])) and float(row["zerr"]) > 0
            and np.isfinite(float(row["z"]))
        )
        result.append({
            "desi_targetid": int(row["targetid"]),
            "desi_ra": float(row["mean_fiber_ra"]), "desi_dec": float(row["mean_fiber_dec"]),
            "desi_z": float(row["z"]), "desi_zerr": float(row["zerr"]),
            "desi_quality_gate": bool(valid_quality),
            "euclid_match": preferred,
            "neighbors": [{
                "id": ident, "ra": mer[ident][0], "dec": mer[ident][1],
                "sep_arcsec": float(sep), "tiles": sorted(tiles[ident]),
            } for sep, ident in neighbors],
        })
    return result


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("desi", type=Path)
    parser.add_argument("mer_manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.mer_manifest.read_bytes())
    with args.desi.open(encoding="utf-8", newline="") as stream:
        ground = list(csv.DictReader(stream))
    if len({r["targetid"] for r in ground}) != len(ground):
        raise ValueError("Expected unique DESI target IDs")
    left, right, bottom, top = manifest["bounds"]
    ground = [row for row in ground if left <= float(row["mean_fiber_ra"]) < right
              and bottom <= float(row["mean_fiber_dec"]) < top]
    mer, tiles = read_mer(manifest["entries"], args.mer_manifest.parent, manifest["bounds"])
    if not mer or not ground:
        raise ValueError("Empty in-bounds source catalog")
    result = associate(ground, mer, tiles)
    output = {"bounds": manifest["bounds"], "mer_unique": len(mer),
              "desi_rows": len(result), "matches": result}
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(output, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
