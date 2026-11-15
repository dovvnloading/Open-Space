"""Search unlabeled GMN trajectories; all thresholds are pairwise."""
import argparse
import csv
import hashlib
from pathlib import Path
from urllib.request import urlopen
import shutil

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

URL = "https://globalmeteornetwork.org/data/traj_summary_data/monthly/traj_summary_monthly_{}.txt"
MONTHS = [f"{y}{m:02}" for y in (2025, 2026) for m in range(1, 13)
          if (y, m) <= (2026, 10)]
COLS = {0: "id", 1: "jd", 2: "utc", 4: "shower", 7: "ra", 9: "dec", 15: "vg",
        63: "lat", 65: "lon", 67: "h", 69: "lat2", 71: "lon2", 73: "h2",
        80: "qc", 84: "nst"}


def unit(lon, lat):
    a, b = np.deg2rad(lon), np.deg2rad(lat)
    return np.column_stack((np.cos(a)*np.cos(b), np.sin(a)*np.cos(b), np.sin(b)))


def read_month(path):
    frame = pd.read_csv(path, sep=";", comment="#", header=None, usecols=list(COLS),
                        skipinitialspace=True, index_col=False).rename(columns=COLS)
    for col in ("id", "utc", "shower"):
        frame[col] = frame[col].str.strip()
    return frame


def xyz(frame):
    # Approximate local Earth as a sphere. Used only for broad separation cuts.
    start = unit(frame.lon, frame.lat) * (6371.0 + frame.h.to_numpy()[:, None])
    end = unit(frame.lon2, frame.lat2) * (6371.0 + frame.h2.to_numpy()[:, None])
    delta = end - start
    length = np.linalg.norm(delta, axis=1)
    return (start + end)/2, delta/length[:, None]


def distinct_paths(frame):
    # Multiple fits of one luminous track can have different starting heights.
    # Collapse connected matches in time, radiant, speed and cross-track distance.
    d = frame.reset_index(drop=True)
    t = pd.to_datetime(d.utc).to_numpy(dtype="datetime64[us]").astype("int64") / 1e6
    r = unit(d.ra, d.dec)
    pos, dr = xyz(d)
    speeds = d.vg.to_numpy()
    parent = np.arange(len(d))

    def root(i):
        while i != parent[i]:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for lag in range(1, len(d)):
        i = np.flatnonzero(t[lag:]-t[:-lag] <= 10)
        if not len(i):
            break
        j = i+lag
        off = pos[i]-pos[j]
        ok = ((np.linalg.norm(off, axis=1) < 80)
              & (np.linalg.norm(np.cross(off, dr[i]), axis=1) < 4)
              & (np.linalg.norm(np.cross(off, dr[j]), axis=1) < 4)
              & (np.linalg.norm(r[i]-r[j], axis=1) < 2*np.sin(np.deg2rad(5)/2))
              & (np.abs(speeds[i]-speeds[j]) < 3))
        for a, b in zip(i[ok], j[ok]):
            parent[root(b)] = root(a)
    d["component"] = [root(i) for i in range(len(d))]
    return (d.sort_values(["nst", "qc", "id"], ascending=[False, False, True])
            .drop_duplicates("component").sort_values("utc").reset_index(drop=True))


def scan(frame, rng, trials, block, radius, regional=False):
    selected = frame.shower.eq("...") & frame.qc.ge(5)
    d = frame[selected].dropna(subset=list(COLS.values())).copy()
    d = d[(d.h > d.h2) & (d.vg > 0)].sort_values("utc").reset_index(drop=True)
    if regional:
        d = distinct_paths(d)
    t = pd.to_datetime(d.utc).to_numpy(dtype="datetime64[us]").astype("int64") / 1e6
    # Entire UTC block or its neighbor could generate a <=10 s coincidence.
    cells = (t // block).astype("int64")
    r = unit(d.ra, d.dec)
    chord = 2*np.sin(np.deg2rad(4)/2)
    features = np.column_stack((r/chord, d.vg/2.2, (cells-cells.min())))
    candidate = cKDTree(features).query_pairs(1.0000000001, p=np.inf, output_type="ndarray")
    i, j = candidate.T
    midpoint, direction = xyz(d)
    off = midpoint[i] - midpoint[j]
    cross_dist = np.minimum(np.linalg.norm(np.cross(off, direction[i]), axis=1),
                            np.linalg.norm(np.cross(off, direction[j]), axis=1))
    ground = unit(d.lon, d.lat)
    sep = np.linalg.norm(ground[i] - ground[j], axis=1)*6371
    valid = ((np.linalg.norm(r[i]-r[j], axis=1) <= chord) & (sep <= radius)
             & (cross_dist >= 10) & (np.abs(d.vg.to_numpy()[i]-d.vg.to_numpy()[j]) <= 2.2))
    # Enumerate triangles and 4-cliques once. Every simulated set of timestamps
    # that can produce a hit is represented in this compatibility graph.
    after = {}
    for a, b in candidate[valid]:
        after.setdefault(int(a), set()).add(int(b))
    groups = {3: [], 4: []}
    for a, neighbors in after.items():
        for b in neighbors:
            both = neighbors & after.get(b, set())
            for c in both:
                groups[3].append((a, b, c))
                groups[4].extend((a, b, c, e) for e in both & after.get(c, set()))
    arrays = {k: np.array(v, dtype=int).reshape(-1, k) for k, v in groups.items()}
    hits = []
    for k, arr in arrays.items():
        for group in arr[np.ptp(t[arr], axis=1) <= 10]:
            hits.append([k, *d.id.to_numpy()[group]])
    # Shuffle the actual timestamps within each UTC block. Keep all observed
    # times, seasonal coverage, and each meteor's joint geometry and velocity.
    # Draw only indices used by possible groups, but sample without replacement
    # in each block, including timestamps of eligible non-candidates.
    all_indices = np.unique(np.concatenate([a.ravel() for a in arrays.values()]))
    compact = {k: np.searchsorted(all_indices, a) for k, a in arrays.items()}
    times = np.zeros((trials, len(all_indices)))
    strata = cells
    if regional:
        # Fixed geographic tiles retain regional activity and camera coverage.
        zone = (np.floor((d.lat.to_numpy()+90)/15)*24
                + np.floor((d.lon.to_numpy()+180)/15)).astype(int)
        strata = cells*400 + zone
    sorted_strata = np.argsort(strata, kind="stable")
    sorted_keys = strata[sorted_strata]
    for key in np.unique(strata[all_indices]):
        slots = np.flatnonzero(strata[all_indices] == key)
        pool = t[sorted_strata[np.searchsorted(sorted_keys, key):
                               np.searchsorted(sorted_keys, key, side="right")]]
        # Per-row random ordering of the timestamp pool is uniform over permutations.
        order = rng.random((trials, len(pool))).argsort(axis=1)[:, :len(slots)]
        times[:, slots] = pool[order]
    controls = np.zeros((trials, 2), dtype=int)
    for k, arr in compact.items():
        for chunk in np.array_split(arr, max(1, (len(arr)+1999)//2000)):
            controls[:, k-3] += (np.ptp(times[:, chunk], axis=2) <= 10).sum(axis=1)
    return len(d), hits, controls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--regional", action="store_true", help="deduplicate paths and shuffle by sky coverage on Earth")
    parser.add_argument("--block", type=int, default=300)
    parser.add_argument("--radius", type=float, default=600)
    parser.add_argument("--trials", type=int, default=1000)
    args = parser.parse_args()
    if args.block < 10 or args.trials < 1:
        parser.error("block must be at least ten seconds and trials must be positive")
    args.cache.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    with (Path(__file__).resolve().parent/"primary"/"months.csv").open() as f:
        frozen = {row["month"]: row["sha256"] for row in csv.DictReader(f)}
    rng = np.random.default_rng(20260824)
    total = np.zeros((args.trials, 2), dtype=int)
    months, hits = [], []
    for month in MONTHS:
        path = args.cache / (month + ".txt")
        if not path.exists() and args.fetch:
            tmp = path.with_suffix(".part")
            with urlopen(URL.format(month), timeout=180) as src, tmp.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            tmp.replace(path)
        with path.open("rb") as f:
            digest = hashlib.file_digest(f, "sha256").hexdigest()
        if digest != frozen[month]:
            print(f"Archive revision in {month}: frozen {frozen[month]}, retrieved {digest}", flush=True)
        frame = read_month(path)
        n, matched, draws = scan(frame, rng, args.trials, args.block, args.radius, args.regional)
        total += draws
        hits.extend(matched)
        months.append(dict(month=month, sha256=digest, rows=len(frame), eligible=n,
                           triples=sum(x[0]==3 for x in matched), quads=sum(x[0]==4 for x in matched)))
        print(months[-1], flush=True)
    with (args.output/"months.csv").open("w") as out:
        writer = csv.DictWriter(out, fieldnames=list(months[0]))
        writer.writeheader()
        writer.writerows(months)
    with (args.output/"groups.csv").open("w") as out:
        writer = csv.writer(out)
        writer.writerow(["multiplicity", "id1", "id2", "id3", "id4"])
        writer.writerows(sorted([*row, *([""] if row[0] == 3 else [])] for row in hits))
    np.savetxt(args.output/"shuffles.csv", total, fmt="%d", delimiter=",",
               header="triples,quadruples", comments="")


if __name__ == "__main__":
    main()
