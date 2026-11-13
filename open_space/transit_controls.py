"""Measure conditional transit recovery across a frozen selection of real TESS backgrounds."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import hashlib
import json
from pathlib import Path

from .tess import injection_experiment
from .tess_archive import atomic_bytes


def run(plan_files, cache, backgrounds=8, workers=2):
    targets = {}
    plans = []
    for path in plan_files:
        raw = path.read_bytes()
        plans.append(hashlib.sha256(raw).hexdigest())
        for target in json.loads(raw)["targets"]:
            if target["role"] == "exploration":
                existing = targets.setdefault(target["tic_id"], target)
                if existing != target:
                    raise ValueError("Conflicting products for a calibration target")
    if not 1 <= backgrounds <= len(targets) or not 1 <= workers <= 32:
        raise ValueError("Invalid number of workers or real stellar backgrounds")
    seed = "open-space-v1-injection-background"
    selection = sorted(
        targets, key=lambda tic: hashlib.sha256(f"{seed}:{tic}".encode()).hexdigest()
    )[:backgrounds]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        trials = list(
            pool.map(
                partial(injection_experiment, cache=cache),
                [targets[tic] for tic in selection],
            )
        )
    by_depth = {}
    for case in trials:
        case["products"] = targets[case["background_tic"]]["products"]
        for trial in case["trials"]:
            depth = trial["injected_depth"]
            row = by_depth.setdefault(
                depth, {"depth": depth, "recovered": 0, "trials": 0}
            )
            row["trials"] += 1
            row["recovered"] += int(trial["recovered"])
    return {
        "experiment": "spoc-transit-injection-backgrounds-v1",
        "plan_hashes": plans,
        "selection": "ascending SHA-256 of seed:TIC among all exploratory targets",
        "seed": seed,
        "backgrounds": backgrounds,
        "available_targets": len(targets),
        "pipeline_stage": "multiplicative dimming in archived corrected photometry before local detrending",
        "depth_summary": list(by_depth.values()),
        "results": trials,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan", action="append", required=True, type=Path)
    p.add_argument("--cache", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--backgrounds", type=int, default=8)
    p.add_argument("--workers", type=int, default=2)
    a = p.parse_args()
    packet = run(a.plan, a.cache, a.backgrounds, a.workers)
    atomic_bytes(
        a.output, (json.dumps(packet, indent=2, allow_nan=False) + "\n").encode()
    )


if __name__ == "__main__":
    main()
