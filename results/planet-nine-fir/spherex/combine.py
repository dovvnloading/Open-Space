"""Combine SPHEREx detector tests using common trial and control identifiers."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def read_json(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    root = parser.parse_args().root
    detectors = [root] + [root / f"D{i}" for i in range(2, 7)]
    trials = []
    eligible_grid = set()
    eligible_controls = set()
    passed_controls = set()
    passed_grid = set()
    total_frames = 0
    individual = []
    for i, directory in enumerate(detectors, 1):
        summary = read_json(directory / "summary.json")
        with gzip.open(directory / "injection_trials.csv.gz", "rt") as f:
            trial = list(csv.DictReader(f))
        trials.append(trial)
        total_frames += summary["exposures"]
        eligible_grid.update(read_json(directory / "eligible_grid_indices.json"))
        eligible_controls.update(read_json(directory / "eligible_control_indices.json"))
        passed_controls.update(read_json(directory / "triggered_control_indices.json"))
        if (directory / "triggers.csv").exists():
            with (directory / "triggers.csv").open() as f:
                passed_grid.update(int(row["grid_index"]) for row in csv.DictReader(f))
        with (directory / "injections.csv").open() as f:
            mag_rows = list(csv.DictReader(f))
        qualifying = [float(row["magnitude_AB"]) for row in mag_rows
                      if int(row["eligible"]) and int(row["recovered"]) / int(row["eligible"]) >= 0.9]
        individual.append(dict(detector=f"D{i}", exposures=summary["exposures"],
                               grid_valid=summary["grid_valid"], grid_triggers=summary["grid_passes"],
                               control_valid=summary["valid_controls"], control_triggers=summary["control_passes"],
                               flux_trials_eligible=summary["eligible_injected_tracks"],
                               faintest_tested_mag_AB_recovery_at_least_90pc=max(qualifying) if qualifying else None))
    ids = [row["grid_index"] for row in trials[0]]
    assert all([row["grid_index"] for row in trial] == ids for trial in trials)
    names = list(trials[0][0])
    assert all(list(trial[0]) == names for trial in trials)
    joint = {}
    for key in names[1:]:
        matrix = np.array([[int(row[key]) for row in trial] for trial in trials])
        joint[key] = int(np.any(matrix, axis=0).sum())
    family_count = read_json(root / "summary.json")["grid_size"]
    coverage = len(eligible_grid) / family_count
    result = dict(exposures=total_frames, grid_size=family_count, grid_valid=len(eligible_grid),
                  grid_passes=len(passed_grid), controls=read_json(root/"summary.json")["controls"],
                  valid_controls=len(eligible_controls), control_passes=len(passed_controls),
                  injected_tracks=len(ids), eligible_injected_tracks=joint.pop("eligible"),
                  detectors=individual, recovered_joint=joint)
    (root/"combined.json").write_text(json.dumps(result, indent=2)+"\n")
    with (root/"combined_injections.csv").open("w", newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["injection", "injected", "eligible", "recovered", "conditional_recovery",
                         "raw_sample_recovery", "sample_recovery_scaled_to_known_coverage"])
        eligible=result["eligible_injected_tracks"]
        for name, n in joint.items():
            writer.writerow([name,len(ids),eligible,n,n/eligible,n/len(ids),coverage*n/eligible])


if __name__ == "__main__":
    main()
