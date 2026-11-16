"""Repeat the holdout rank comparison on equally exposed H1 epochs at every shift.

Uses the released, unrounded sparse tail rows and per-detector searchable masks.
This post-design diagnostic does not change the search statistic or rate limit.
"""
import argparse
import csv
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np


def aligned_masks(keep_h, keep_l, shifts):
    n = len(keep_h)
    common = keep_h.copy()
    for shift in shifts:
        lo, hi = max(0, -shift), min(n, n-shift)
        common[:lo] = False
        common[hi:] = False
        common[lo:hi] &= keep_l[lo+shift:hi+shift]
    return common


def calculate(directory):
    campaign = json.loads((directory / "campaign.json").read_text())
    shifts = [r[0] for r in campaign["campaign_slides"]]
    maxima = dict.fromkeys(shifts)
    epochs = []
    for shard in campaign["monthly_ledgers"]:
        with ZipFile(directory / shard["file"]) as z:
            for member in z.namelist():
                if not member.endswith("/epoch.json"):
                    continue
                data = json.loads(z.read(member))
                masks = []
                for detector in ("H1", "L1"):
                    keep = np.zeros(4096, dtype=bool)
                    for lo, hi in data["masks"][detector]["searchable_runs"]:
                        keep[lo:hi] = True
                    masks.append(keep)
                common = aligned_masks(*masks, shifts)
                epochs.append({"gps":data["gps"], "common_H_seconds":int(common.sum())})
                tail = csv.DictReader(io.StringIO(z.read(member.replace("epoch.json","tail.csv")).decode()))
                for row in tail:
                    t = float(row["H1_peak_gps"])-data["gps"]
                    s, score = int(row["shift_s"]), float(row["score"])
                    if common[int(t)] and (maxima[s] is None or score > maxima[s]):
                        maxima[s] = score
    if len(epochs) != campaign["n_frames"] or maxima[0] is None:
        raise ValueError("Missing epochs or foreground below sparse-tail floor")
    n = sum(x is not None and x >= maxima[0] for s, x in maxima.items() if s)
    return {"common_H_seconds_per_campaign":sum(e["common_H_seconds"] for e in epochs),
            "foreground_maximum":maxima[0], "shifted_campaigns_ge_foreground":n,
            "empirical_plus_one_fraction":(n+1)/len(shifts),
            "maxima_by_shift":maxima, "epochs":epochs, "tail_storage_floor":3.0,
            "description":"All 129 alignments have precisely the same eligible H1 seconds. Each alignment uses a different L1 epoch. Null per-shift maxima indicate scores below the 3.0 storage floor."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=Path(__file__).parent/"results"/"holdout")
    args = parser.parse_args()
    output = calculate(args.directory)
    (args.directory/"common_exposure.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    print({k:v for k,v in output.items() if k not in ("maxima_by_shift","epochs")})
