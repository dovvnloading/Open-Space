"""Rebuild the ATNF and two Pushchino screens, including the field-pulsar controls."""
import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import tarfile
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
CATALOG = "https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.v2.9.0.tar.gz"
DB_HASH = "bac32471fef3e60cfd734f1db85503a9bd0a565018cb98edff54475f87eb642f"
HARMONICS = np.array([1, 1.5, 2, 3, 4, 5, 8])


def sexagesimal(text):
    sign = -1 if text.startswith("-") else 1
    return sign * sum(float(v) / 60**i for i, v in enumerate(text.lstrip("+-").split(":")))


def load_database(cache):
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / "psrcat-2.9.0.db"
    if not path.exists():
        with urlopen(CATALOG, timeout=60) as response:
            archive = response.read()
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as tar:
            data = tar.extractfile("psrcat_tar/psrcat.db").read()
    else:
        data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != DB_HASH:
        raise ValueError("Unexpected ATNF database digest")
    if not path.exists():
        path.write_bytes(data)
    objects = []
    for block in data.decode("latin1").split("\n@"):
        row = {}
        for line in block.splitlines():
            fields = line.split()
            if len(fields) >= 2 and not line.startswith(("#", "-")):
                row.setdefault(fields[0], fields[1])
        if "PSRJ" in row:
            objects.append(row)
    return objects


def coordinates(r):
    if "RAJ" in r and "DECJ" in r:
        return 15 * sexagesimal(r["RAJ"]), sexagesimal(r["DECJ"])
    lon, lat, eps = np.deg2rad([float(r["ELONG"]), float(r["ELAT"]), 23.4393])
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon) * np.cos(eps) - np.sin(lat) * np.sin(eps)
    z = np.cos(lat) * np.sin(lon) * np.sin(eps) + np.sin(lat) * np.cos(eps)
    return np.rad2deg(np.arctan2(y, x)) % 360, np.rad2deg(np.arcsin(z))


def write_rows(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def atnf_screen(objects, out):
    rows = [r for r in objects if "DM" in r and ("P0" in r or "F0" in r)
            and (("RAJ" in r and "DECJ" in r) or ("ELONG" in r and "ELAT" in r))]
    period = np.array([float(r["P0"]) if "P0" in r else 1 / float(r["F0"]) for r in rows])
    dm = np.array([float(r["DM"]) for r in rows])
    ra, dec = np.deg2rad([coordinates(r) for r in rows]).T
    sky = np.column_stack((np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)))
    i, j = cKDTree(sky).query_pairs(2*math.sin(math.radians(1.1)/2), output_type="ndarray").T
    sep = np.rad2deg(2*np.arcsin(np.linalg.norm(sky[i]-sky[j], axis=1)/2))
    close = (np.abs(dm[i] - dm[j]) <= 3)
    globular = np.array(["GC:" in r.get("ASSOC", "") for r in rows])
    subset = close & ~globular[i] & ~globular[j]
    ratio = np.maximum(period[i],period[j]) / np.minimum(period[i],period[j])
    errors = np.abs(ratio[:,None] / HARMONICS - 1)
    which = np.argmin(errors, axis=1)
    matches = []
    for k in np.flatnonzero(close & (errors.min(axis=1) <= 5e-4)):
        a, b = i[k], j[k]
        matches.append(dict(a=rows[a]["PSRJ"], b=rows[b]["PSRJ"], separation_deg=f"{sep[k]:.6f}",
                            p_a_s=f"{period[a]:.12g}", p_b_s=f"{period[b]:.12g}",
                            dm_a=dm[a], dm_b=dm[b], ratio=HARMONICS[which[k]],
                            residual_ppm=f"{1e6*errors[k,which[k]]:.3f}",
                            gc_flag=bool(globular[a] or globular[b])))
    write_rows(out / "screen.csv", sorted(matches, key=lambda m:(m["a"],m["b"])))

    # Shuffle only field periods, locally in declination and log DM. Retain each real close
    # pair, common sky/DM environment, total trials, and empirical period distribution.
    groups = defaultdict(list)
    for k in np.flatnonzero(~globular):
        groups[(math.floor(math.degrees(dec[k])/15),math.floor(math.log10(max(1,dm[k]))/0.25))].append(k)
    groups = [np.array(g) for g in groups.values() if len(g)>1]
    left, right = i[subset], j[subset]
    rng = np.random.default_rng(20261115)
    null = []
    for _ in range(1000):
        sim = period.copy()
        for g in groups:
            sim[g] = rng.permutation(period[g])
        r = np.maximum(sim[left],sim[right])/np.minimum(sim[left],sim[right])
        null.append(int(np.count_nonzero(np.min(np.abs(r[:,None]/HARMONICS - 1),axis=1)<=5e-4)))
    field = np.flatnonzero(~globular)
    tree = cKDTree(sky[field])
    shifted_counts = {}
    # Shift every source, then compare it with the unchanged field catalog. These are
    # directed cross-catalog matches; count both directions in the zero-shift baseline.
    for shift in (-12, -8, -4, 4, 8, 12):
        angle = np.deg2rad(shift)
        shifted = sky[field].copy()
        shifted[:,0] = np.cos(dec[field]) * np.cos(ra[field] + angle)
        shifted[:,1] = np.cos(dec[field]) * np.sin(ra[field] + angle)
        neighborhoods = tree.query_ball_point(shifted, 2*math.sin(math.radians(1.1)/2))
        shifted_hits = 0
        for source, neighbors in zip(field, neighborhoods):
            candidates = field[neighbors]
            keep = (candidates!=source) & (np.abs(dm[candidates]-dm[source])<=3)
            others = candidates[keep]
            rr = np.maximum(period[others],period[source])/np.minimum(period[others],period[source])
            shifted_hits += int(np.count_nonzero(np.min(np.abs(rr[:,None]/HARMONICS - 1),axis=1)<=5e-4))
        shifted_counts[str(shift)] = shifted_hits
    stats = dict(catalog=len(objects), eligible=len(rows), geometric_pairs=len(i),
                 geometric_dm_pairs=int(close.sum()), field_geometric_dm_pairs=len(left),
                 observed_field_pairs=sum(not m["gc_flag"] for m in matches),
                 null_trials=len(null), null_min=min(null),null_median=float(np.median(null)),
                 null_max=max(null),null_mean=float(np.mean(null)),
                 shifted_ra_directed_counts=shifted_counts)
    return stats


def pushchino_screen(out):
    with (HERE/"pushchino-2020.csv").open() as f, (HERE/"pushchino-2024.csv").open() as g:
        old,new = list(csv.DictReader(f)),list(csv.DictReader(g))
    joined = []
    for a,b in itertools.product(old,new):
        # Published coordinate errors are boxes. RA units are seconds of sidereal time.
        dra = (sexagesimal(a["ra_hms"])-sexagesimal(b["ra_hms"]))*3600
        dra = (dra+43200)%86400-43200
        ddec = (sexagesimal(a["dec_dm"])-sexagesimal(b["dec_dm"]))*60
        checks = [(dra,"ra_error_s_time"), (ddec,"dec_error_arcmin"),
                  (float(a["period_s"])-float(b["period_s"]),"period_error_s"),
                  (float(a["dm"])-float(b["dm"]),"dm_error")]
        passes = [abs(d)<=float(a[col])+float(b[col]) for d,col in checks]
        joined.append(dict(older=a["name"],later=b["name"],ra_delta_s=f"{dra:.3f}",
                           dec_delta_arcmin=f"{ddec:.3f}",period_delta_s=f"{checks[2][0]:.7f}",
                           dm_delta=f"{checks[3][0]:.3f}",boxes_overlap=all(passes[:2]),
                           all_intervals_overlap=all(passes)))
    write_rows(out/"pushchino-pairs.csv",joined)
    return dict(pairs=len(joined),sky_overlaps=sum(i["boxes_overlap"] for i in joined),
                four_way_overlaps=sum(i["all_intervals_overlap"] for i in joined))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    objects = load_database(args.cache)
    summary = {"atnf":atnf_screen(objects,args.out),"pushchino":pushchino_screen(args.out)}
    (args.out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary))
