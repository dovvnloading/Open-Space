"""Shuffle whole-spectrum redshift labels to measure rest-frame recurrence under reassignment."""

import argparse
from collections import Counter
import hashlib
from itertools import permutations, product
import json
import math
from pathlib import Path
import random

from .dataio import canonical, digest
from .recurrence import _clusters, packet_files, summarize
from .spectra import CATALOG_VERSION, _identify


def _statistics(groups):
    counts = [g["target_count"] for g in groups]
    return {
        "largest_group": max(counts, default=0),
        "groups_with_at_least_four_targets": sum(n >= 4 for n in counts),
        "groups_with_at_least_two_targets": len(counts),
    }


def _records(packets):
    records = []
    for packet in sorted(packets, key=lambda p: p["observation"]["observation_id"]):
        obs, scan = packet["observation"], packet["stages"]["detection"]
        z = obs["redshift"]
        if not math.isfinite(z) or z <= -1 or z != scan["redshift"]:
            raise ValueError("Invalid or inconsistent source redshift")
        observed = []
        for index, f in enumerate(scan["features"], 1):
            if (not all(math.isfinite(f[k]) and f[k] > 0 for k in
                        ("wavelength", "score", "pixel_spacing_angstrom")) or
                    f["direction"] not in ("absorption", "emission")):
                raise ValueError("Invalid detected feature")
            identity = _identify(f, z)
            if (identity["classification"] != f["classification"] or
                    not math.isclose(identity["rest_wavelength"], f["rest_wavelength"], rel_tol=1e-12)):
                raise ValueError("Identity redshift does not reproduce archived feature classification or wavelength")
            observed.append({k: f[k] for k in ("wavelength", "direction", "score", "pixel_spacing_angstrom",
                                                "classification")} | {"feature_id": f"F{index:04d}"})
        records.append((obs, observed))
    return records


def _reassign(records, assigned, target, donor):
    obs, features = records[target]
    z = assigned[donor]
    hits, frozen_hits = [], []
    for feature in features:
        classified = _identify(feature, z)
        hit = {k: classified[k] for k in ("wavelength", "rest_wavelength", "direction", "score", "feature_id")}
        hit.update(target=target, plate=obs["plate"], mjd=obs["mjd"],
                   source_sha256=obs["source_sha256"], observation_id=obs["observation_id"])
        if classified["classification"] == "unidentified":
            hits.append(hit)
        if feature["classification"] == "unidentified":
            frozen_hits.append(hit)
    return hits, frozen_hits


def _assignments(blocks, draws, seed):
    """Uniformly permute donor IDs within each block, keeping spectra intact.

    Small spaces are enumerated once including identity. Otherwise draw
    independent uniform assignments with replacement, allowing fixed points.
    """
    size = math.prod(math.factorial(len(b)) for b in blocks)
    n = sum(len(b) for b in blocks)
    if size <= draws:
        # Each block permutation is short when the *combined* space fits the budget.
        combinations = product(*(permutations(b) for b in blocks))
    else:
        rng = random.Random(seed)
        combinations = (tuple(rng.sample(b, len(b)) for b in blocks) for _ in range(draws))

    def generate():
        for assignment in combinations:
            order = list(range(n))
            for block, donors in zip(blocks, assignment):
                for target, donor in zip(block, donors):
                    order[target] = donor
            yield order

    return size, generate()


def _evaluate(records, remapped, blocks, draws, seed, tolerance, observed, fixed_sensitivity):
    modes = ["reapply_catalog", "freeze_original_shortlist"] if fixed_sensitivity else ["reapply_catalog"]
    results = {mode: {k: Counter() for k in observed} for mode in modes}
    identity = list(range(len(records)))
    sequence = hashlib.sha256()
    size, assignments = _assignments(blocks, draws, seed)
    count = 0
    for order in assignments:
        sequence.update((canonical(order) + "\n").encode())
        for mode_index, mode in enumerate(modes):
            hits = [hit for target, donor in enumerate(order) for hit in remapped[target][donor][mode_index]]
            groups = _clusters(hits, "rest_wavelength", tolerance, 2)
            stats = _statistics(groups) | {"features_after_catalog_screen": len(hits)}
            if order == identity and observed != stats:
                raise ValueError("Identity reassignment changed the original recurrence summary")
            for name, value in stats.items():
                results[mode][name][value] += 1
        count += 1
    exact = size <= draws
    distributions = {}
    for mode, histograms in results.items():
        distributions[mode] = {}
        for name, hist in histograms.items():
            entry = {"histogram": {str(k): hist[k] for k in sorted(hist)}, "observed": observed[name]}
            if name in ("largest_group", "groups_with_at_least_four_targets"):
                ge = sum(num for value, num in hist.items() if value >= observed[name])
                entry["draws_at_least_observed"] = ge
                entry["upper_tail_fraction"] = ge / count if exact else (ge + 1) / (count + 1)
            distributions[mode][name] = entry
    return {
        "mode": "exact" if exact else "monte_carlo_with_replacement", "assignments": count,
        "possible_index_assignments": size, "assignment_sequence_sha256": sequence.hexdigest(),
        "blocks": [[records[i][0]["observation_id"] for i in b] for b in blocks],
        "distributions": distributions,
    }


def experiment(packets, draws=9999, seed=6497, tolerance_kms=150.0):
    """Compare observed grouping with spectrum-wise randomized redshift assignments."""
    if isinstance(draws, bool) or not isinstance(draws, int) or draws < 1:
        raise ValueError("draws must be a positive integer")
    packets = list(packets)
    baseline = summarize(packets, tolerance_kms)
    if len(packets) < 2 or baseline["distinct_targets"] != len(packets):
        raise ValueError("Require at least two distinct sky targets, each with exactly one completed packet")
    if baseline["compatibility"]["catalog"] != CATALOG_VERSION:
        raise ValueError("Catalog replay requires the same version used in archived scans")
    records = _records(packets)
    z = [obs["redshift"] for obs, _ in records]
    if len(set(z)) < 2:
        raise ValueError("Require at least two distinct redshifts")
    if len({(obs["plate"], obs["mjd"]) for obs, _ in records}) != 1:
        raise ValueError("Supply one plate/MJD visit per experiment")
    observed = _statistics(baseline["rest_frame_groups"]) | {
        "features_after_catalog_screen": baseline["unidentified_features"]}
    indices = list(range(len(records)))
    # Redshift and source define a small table for this single-visit experiment.
    # Reusing the table avoids repeating catalog lookup during each draw.
    remapped = [[_reassign(records, z, i, j) for j in indices] for i in indices]
    identity_hits = [hit for i in indices for hit in remapped[i][i][0]]
    if _clusters(identity_hits, "rest_wavelength", tolerance_kms, 2) != baseline["rest_frame_groups"]:
        raise ValueError("Identity assignment disagrees with original recurrence groups")
    ordered = sorted(indices, key=lambda i: (z[i], i))
    nearby = [ordered[i:i + 3] for i in range(0, len(ordered), 3)]
    inputs = [{"observation_id": obs["observation_id"], "source_sha256": obs["source_sha256"],
               "csv_sha256": obs["csv_sha256"], "redshift": obs["redshift"],
               "detected_features": len(features), "feature_vector_sha256": digest(canonical(features).encode())}
              for obs, features in records]
    return {
        "schema_version": 1, "protocol": "whole-spectrum-redshift-reassignment-v1",
        "claim_status": "conditional_control", "draw_budget": draws, "seed": seed,
        "tolerance_kms": tolerance_kms, "compatibility": baseline["compatibility"],
        "implementation_sha256": {name: digest(Path(__file__).with_name(name).read_bytes())
                                  for name in ("redshift_null.py", "recurrence.py", "spectra.py", "dataio.py")},
        "inputs": inputs, "visit": {"plate": records[0][0]["plate"], "mjd": records[0][0]["mjd"]},
        "observed": observed, "observer_frame_original_shortlist": _statistics(baseline["observer_frame_control_groups"]),
        "unrestricted": _evaluate(records, remapped, [indices], draws, seed, tolerance_kms, observed, True),
        "neighboring_redshift_triples": _evaluate(records, remapped, nearby, draws, seed, tolerance_kms, observed, False),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="Completed evidence packets, reports directory or campaign")
    parser.add_argument("--draws", type=int, default=9999, help="Maximum assignments per experiment")
    parser.add_argument("--seed", type=int, default=6497)
    parser.add_argument("--tolerance-kms", type=float, default=150.0)
    args = parser.parse_args(argv)
    try:
        packets = [json.loads(p.read_bytes()) for p in packet_files(args.inputs)]
        result = experiment(packets, args.draws, args.seed, args.tolerance_kms)
        print(json.dumps(result, indent=2, allow_nan=False))
    except (ValueError, TypeError, KeyError, OSError) as exc:
        parser.exit(1, f"redshift-null failed: {exc}\n")


if __name__ == "__main__":
    main()
