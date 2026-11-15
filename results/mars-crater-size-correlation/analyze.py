"""Test spatial correlation of measured diameters within fresh Martian crater clusters.

Run with Python 3, NumPy, SciPy and openpyxl. See README.md for source and usage.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from openpyxl import load_workbook
from scipy.spatial.distance import cdist, pdist
from scipy.stats import binomtest, rankdata

SOURCE_SHA256 = "f32fd7ece38d8ebb9a8874651d14f09c450278f8e58a27858c12f88219651aed"
RADIUS_M = 3389500.0


def percentile_ranks(values):
    return (rankdata(values) - 1) / (len(values) - 1)


def equal_count(values, groups):
    order = np.argsort(values, kind="stable")
    labels = np.empty(len(values), dtype=int)
    labels[order] = np.arange(len(values)) * groups // len(values)
    return labels


def expected_edge_gap(marks, labels, start, end):
    """Exact mean under exchangeability without replacement inside each stratum."""
    levels = np.unique(labels)
    gaps = np.empty((max(levels) + 1, max(levels) + 1))
    for first in levels:
        left = marks[labels == first]
        for second in levels:
            right = marks[labels == second]
            pairs = len(left) * len(right) - (len(left) if first == second else 0)
            gaps[first, second] = (np.abs(left[:, None] - right).sum() / pairs
                                    if pairs else np.nan)
    return float(gaps[labels[start], labels[end]].mean())


def site_statistics(xy_m, diameter_m, null_draws=None, rng=None):
    """One site's equal-weight summaries. Caller applies the diameter threshold."""
    marks = percentile_ranks(diameter_m)
    n = len(marks)
    if n < 12:
        raise ValueError("At least 12 craters required.")
    xy = np.asarray(xy_m, dtype=float)
    centered = xy - xy.mean(axis=0)
    axes = centered @ np.linalg.svd(centered, full_matrices=False)[2].T
    # Set the arbitrary sign consistently, including for unequal-sized bins.
    for column in range(2):
        farthest = np.argmax(np.abs(axes[:, column]))
        if axes[farthest, column] < 0:
            axes[:, column] *= -1
    groups = min(4, n // 4)
    along = equal_count(axes[:, 0], groups)
    distances = cdist(xy, xy)
    np.fill_diagonal(distances, np.inf)
    neighbors = np.argsort(distances, axis=1, kind="stable")
    start = np.arange(n)
    end = neighbors[:, 0]
    if np.any(distances[start, end] == 0):
        raise ValueError("Coincident crater centers.")
    if null_draws is not None:
        # Independent within-field, within-bin permutations, summed by field.
        # Batch to bound memory at the largest observed site (over 2,000 craters).
        for offset in range(0, len(null_draws), 64):
            size = min(64, len(null_draws) - offset)
            shuffled = np.empty((size, n))
            for label in range(groups):
                cell = np.flatnonzero(along == label)
                keys = rng.random((size, len(cell)))
                shuffled[:, cell] = marks[cell[np.argsort(keys, axis=1)]]
            null_draws[offset:offset + size] += np.abs(shuffled - shuffled[:, end]).mean(axis=1)
    stats = {}

    def add(name, labels, edges, other=end):
        if not edges.any():
            return
        u, v = start[edges], other[edges]
        observed = float(np.abs(marks[u] - marks[v]).mean())
        expected = expected_edge_gap(marks, labels, u, v)
        stats[name] = (len(u), observed, expected)

    all_edges = np.ones(n, dtype=bool)
    add("global", np.zeros(n, dtype=int), all_edges)
    add("along", along, all_edges)
    # Small spatial cells and local crowding are separate sensitivity tests.
    # Bisect each along-track block at its cross-track or local-density median.
    for name, variable in [
        ("along_cross", axes[:, 1]),
        ("along_density", distances[start, neighbors[:, 3]]),
    ]:
        labels = along * 2
        for group in range(groups):
            cell = along == group
            labels[cell] += equal_count(variable[cell], 2)
        add(name, labels, all_edges)
    # These masks use positions and the sitewide maximum only. Reassigning
    # marks therefore cannot alter which edges survive.
    clearance = distances[start, end] / np.max(diameter_m)
    for scale in (2, 4):
        use = clearance >= scale
        if use.sum() >= 6:
            add(f"separated_{scale}dmax", along, use)
    isolated = clearance > 1
    edges = isolated & isolated[end]
    if edges.sum() >= 6:
        # Reassign sizes only among isolated positions; every reassignment is
        # physically disjoint from every other mapped rim. Keep original bins.
        restricted = along.copy()
        restricted[~isolated] += groups
        add("isolated_exchange", restricted, edges)
    for rank in (2, 4, 8):
        add(f"neighbor_{rank}", along, all_edges, neighbors[:, rank - 1])
    return stats


def cluster_rows(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("Input differs from Neidhart et al. record 5760176.")
    with zipfile.ZipFile(path) as archive:
        data = archive.read("cluster_data/crater_clusters.csv").decode("latin1")
        for field in csv.DictReader(io.StringIO(data)):
            image = field["HiRiseID"]
            name = f"cluster_data/data_set/{image}.xlsx"
            sheet = load_workbook(io.BytesIO(archive.read(name)), read_only=True,
                                  data_only=True).active
            rows = iter(sheet.values)
            header = next(rows)
            indices = [header.index(key) for key in ("Diam_km", "x_coord", "y_coord")]
            values = np.array([[r[c] for c in indices] for r in rows
                               if all(isinstance(r[c], (int, float)) for c in indices)])
            if len(values) != int(field["N_c"]):
                raise ValueError(f"Count mismatch at {image}")
            diameter_m, lon, lat = values.T
            if not (np.all(diameter_m > 0) and np.all(np.abs(lat) <= 90)):
                raise ValueError(f"Invalid values at {image}")
            # Several workbooks use local longitude conventions. Only angular
            # differences inside each workbook enter the measurement.
            dx = (lon - lon[0] + 180) % 360 - 180
            dy = lat - lat[0]
            xy = np.deg2rad(np.column_stack((dx * np.cos(np.deg2rad(lat.mean())), dy))) * RADIUS_M
            yield field, xy, 1000 * diameter_m


def summarize(rows, rng):
    observed = np.array([r["observed"] for r in rows])
    expected = np.array([r["expected"] for r in rows])
    delta = expected - observed
    boot = rng.choice(delta, (10000, len(rows)), replace=True).mean(axis=1)
    signs = delta != 0
    return dict(fields=len(rows), craters=sum(r["n"] for r in rows),
                directed_edges=sum(r["edges"] for r in rows),
                positive=int((delta > 0).sum()),
                observed=float(observed.mean()), expected=float(expected.mean()),
                difference=float(delta.mean()),
                difference_ci95=np.quantile(boot, [.025, .975]).tolist(),
                sign_test_p_two_sided=float(binomtest(int((delta > 0).sum()),
                                                      int(signs.sum())).pvalue))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    results = []
    science = []
    null_draws = np.zeros(9999)
    null_rng = np.random.default_rng(440321)
    for field, xy, diam in cluster_rows(args.archive):
        image = field["HiRiseID"]
        orbit = int(image.split("_")[1])
        reference = float(field["d_med"])
        separation = float(np.median(pdist(xy)))
        science.append(dict(id=image, n=len(diam), published_median_m=reference,
                            measured_median_m=separation))
        for threshold in (1, 1.5, 2, 3):
            use = diam >= threshold
            if use.sum() < 12:
                continue
            metrics = site_statistics(xy[use], diam[use], null_draws if threshold == 1 else None, null_rng)
            for method, (edges, obs, exp) in metrics.items():
                results.append(dict(id=image, site_id=field["UniqueSiteID"],
                                    orbit=orbit, threshold_m=threshold, n=int(use.sum()),
                                    control=method, edges=edges, observed=obs, expected=exp))
    for name, rows in [("field_statistics.csv", results), ("geometry_audit.csv", science)]:
        with (args.output / name).open("w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    rng = np.random.default_rng(5760176)
    report = dict(fields_in_archive=len(science), craters_in_archive=sum(r["n"] for r in science),
                  source_sha256=SOURCE_SHA256, summaries={})
    for cut in (1, 1.5, 2, 3):
        for mode in sorted({r["control"] for r in results}):
            rows = [r for r in results if r["threshold_m"] == cut and r["control"] == mode]
            if rows:
                report["summaries"][f"{cut}m_{mode}"] = summarize(rows, rng)
    for label, pred in [
        ("orbit_early", lambda r: r["orbit"] < 40000),
        ("orbit_late", lambda r: r["orbit"] >= 40000),
        ("exclude_geometry_outlier", lambda r: r["id"] != "ESP_012047_1755"),
        ("under_50_members", lambda r: r["n"] < 50),
        ("50plus_members", lambda r: r["n"] >= 50),
    ]:
        rows = [r for r in results if r["threshold_m"] == 1 and r["control"] == "along" and pred(r)]
        report["summaries"][label] = summarize(rows, rng)
    primary = report["summaries"]["1m_along"]
    null_draws /= primary["fields"]
    exceed = int(sum(null_draws <= primary["observed"]))
    report["along_track_permutation"] = dict(draws=len(null_draws), seed=440321,
        at_or_below_observed=exceed, p_lower_add_one=(exceed + 1) / (len(null_draws) + 1),
        mean=float(null_draws.mean()), sd=float(null_draws.std(ddof=1)),
        central_95=np.quantile(null_draws, [.025, .975]).tolist())
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
