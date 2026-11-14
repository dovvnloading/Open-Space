"""Replay the sky, dispersion and morphology comparison on frozen CHIME measurements."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
SOURCE = "https://assets.chime-frb.ca/catalogs/catalog2/releases/20251231/data/chimefrbcat2.csv"
FOCUS = ("FRB20200511C", "FRB20211128E")
FIELDS = ("tns_name", "repeater_name", "ra", "dec", "gb", "dm_fitb", "dm_fitb_err",
          "width_fitb", "peak_freq", "low_freq", "high_freq", "mjd_inf", "bonsai_snr",
          "excluded_flag", "sidelobe_flag", "citizen_science_flag")
REQUIRED = ("ra", "dec", "gb", "dm_fitb", "width_fitb", "peak_freq", "low_freq", "high_freq",
            "mjd_inf", "bonsai_snr")


def reduce_catalog(raw):
    """Keep the archive's first component of each event; it lists components in time order."""
    records = {}
    for row in csv.DictReader(io.StringIO(raw.decode("utf-8"))):
        records.setdefault(row["event_id"], row)
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records.values())
    return stream.getvalue().encode()


def sky_vectors(rows):
    ra = np.deg2rad([float(r["ra"]) for r in rows])
    dec = np.deg2rad([float(r["dec"]) for r in rows])
    return np.column_stack((np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)))


def eligible(row):
    return (not row["repeater_name"]
            and all(row[k] == "0" for k in ("excluded_flag", "sidelobe_flag", "citizen_science_flag"))
            and all(np.isfinite(float(row[k])) for k in REQUIRED)
            and float(row["width_fitb"]) > 0)


def analyze(catalog, trials):
    rows = [r for r in catalog if eligible(r)]
    names = np.array([r["tns_name"] for r in rows])
    cols = {k: np.array([float(r[k]) for r in rows]) for k in REQUIRED}
    sky = sky_vectors(rows)
    limit = 2 * np.sin(np.deg2rad(1) / 2)
    pairs = cKDTree(sky).query_pairs(limit, output_type="ndarray")
    i, j = pairs.T
    separation = np.rad2deg(2 * np.arcsin(np.linalg.norm(sky[i] - sky[j], axis=1) / 2))
    pair = np.flatnonzero(((names[i] == FOCUS[0]) & (names[j] == FOCUS[1]))
                         | ((names[j] == FOCUS[0]) & (names[i] == FOCUS[1]))).item()
    angular = 0.25 ** 2 + separation ** 2
    dm = cols["dm_fitb"]
    frequency = cols["peak_freq"]
    log_width = np.log(cols["width_fitb"])

    def score(permutation):
        x, y = permutation[i], permutation[j]
        dispersion = 1 + np.abs(dm[x] - dm[y])
        spectrum = 1 + np.abs(frequency[x] - frequency[y]) / 25
        duration = 1 + np.abs(log_width[x] - log_width[y]) / 0.25
        morphology = spectrum * duration
        return np.array((angular * dispersion * morphology, angular * dispersion,
                         (1 + separation ** 2) * dispersion * morphology))

    identity = np.arange(len(rows))
    original = score(identity)
    primary = original[0]
    score_names = ["position_dm_shape", "position_dm_only", "one_degree_floor_with_shape"]
    c = cols
    dec_bin = np.floor((c["dec"] + 90) / 10).astype(int)
    lat_bin = np.floor(np.abs(c["gb"]) / 20).astype(int)
    time_bin = np.floor((c["mjd_inf"] - 58000) / 500).astype(int)
    tests = {"dec": dec_bin, "dec_latitude": dec_bin * 10 + lat_bin,
             "dec_latitude_epoch": dec_bin * 1000 + lat_bin * 10 + time_bin}
    simulations = {}
    for key, groups in tests.items():
        rng = np.random.default_rng(140726)
        strata = [np.flatnonzero(groups == s) for s in np.unique(groups)]
        hits = np.zeros(len(original), dtype=int)
        for _ in range(trials):
            p = identity.copy()
            for group in strata:
                p[group] = rng.permutation(group)
            hits += score(p).min(axis=1) <= original[:, pair]
        simulations[key] = dict(zip(score_names, hits.tolist()))
    focus_rows = [next(row for row in rows if row["tns_name"] == name) for name in FOCUS]
    a, b = focus_rows
    differences = {key: abs(float(a[key]) - float(b[key])) for key in
                   ("dm_fitb", "peak_freq", "mjd_inf")}
    differences["combined_formal_dm_err"] = np.hypot(float(a["dm_fitb_err"]), float(b["dm_fitb_err"]))
    dm_center = (float(a["dm_fitb"]) + float(b["dm_fitb"])) / 2
    coords = [r for r in catalog if all(np.isfinite(float(r[k])) for k in ("ra", "dec", "dm_fitb"))]
    ref = sky_vectors([a])[0]
    dist = np.rad2deg(2 * np.arcsin(np.minimum(np.linalg.norm(sky_vectors(coords)-ref, axis=1)/2,1)))
    nearby = [r["tns_name"] for p, r in zip(dist, coords)
              if p < 5 and abs(float(r["dm_fitb"]) - dm_center) < 15]
    known_dm = {}
    for row in catalog:
        if row["repeater_name"] and np.isfinite(float(row["dm_fitb"])):
            known_dm.setdefault(row["repeater_name"], []).append(float(row["dm_fitb"]))
    highest = max(known_dm, key=lambda n: np.median(known_dm[n]))
    return {
        "archive_events": len(catalog), "events_searched": len(rows), "pairs_within_one_degree": len(i),
        "focus": list(FOCUS), "separation_arcseconds": float(separation[pair]*3600),
        "absolute_differences": differences, "focus_ranks_in_scores":
        dict(zip(score_names,(original <= original[:,[pair]]).sum(axis=1).tolist())),
        "trials_per_stratification": trials, "catalog_minimum_as_small_as_focus": simulations,
        "top_ten_primary": [[names[i[t]], names[j[t]], round(primary[t], 6)]
                            for t in np.argsort(primary)[:10]],
        "full_catalog_5deg_15DM_check": nearby,
        "largest_median_dm_of_labelled_repeaters": [highest, float(np.median(known_dm[highest]))],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-archive", action="store_true")
    parser.add_argument("--trials", type=int, default=10000)
    args = parser.parse_args()
    if args.trials <= 0:
        parser.error("--trials must be positive")
    receipt = json.loads((ROOT / "sources.json").read_text())
    snapshot = (ROOT / "chime-events.csv").read_bytes()
    if hashlib.sha256(snapshot).hexdigest() != receipt["event_snapshot_sha256"]:
        raise ValueError("Event snapshot does not match its receipt")
    if args.verify_archive:
        with urlopen(SOURCE, timeout=90) as response:
            raw = response.read()
        if hashlib.sha256(raw).hexdigest() != receipt["catalog_sha256"] or reduce_catalog(raw) != snapshot:
            raise ValueError("Published catalog or reduced event table differs from the receipt")
    print(json.dumps(analyze(list(csv.DictReader(io.StringIO(snapshot.decode()))), args.trials), indent=2))
