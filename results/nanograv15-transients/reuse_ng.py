"""Restore unchanged narrowband/wideband reductions on an external-validation replay."""

import json
from pathlib import Path
import shutil
import sys


def restore(dataset, shard):
    if dataset not in ("nb", "wb") or not 0 <= shard < 4:
        raise ValueError("invalid reduction partition")
    previous, current = Path("previous"), Path("reduced")
    receipt = f"coverage-{dataset}-{shard}.json"
    rows = json.loads((previous / receipt).read_text())
    manifest = json.loads(Path("raw/inputs.json").read_text())
    prior = {row["pulsar"]: row for row in rows}
    selected = [row for row in manifest["sources"] if row["dataset"] == dataset][shard::4]
    if set(prior) != {row["pulsar"] for row in selected}:
        raise ValueError("pulsar partition changed")
    for row in selected:
        old = prior[row["pulsar"]]
        for key in ("par", "tim", "par_sha256", "tim_sha256", "include_sha256"):
            if old.get(key) != row.get(key):
                raise ValueError(f"source changed for {dataset}/{row['pulsar']}: {key}")
        if old["status"] == "ok" and not (previous / dataset / (row["pulsar"] + ".npz")).is_file():
            raise ValueError("missing reduced data product")
    shutil.copytree(previous / dataset, current / dataset)
    for name in (receipt, f"runtime-{dataset}-{shard}.json"):
        shutil.copyfile(previous / name, current / name)
    shutil.copyfile(previous / f"software-{shard}.txt", current / f"source-software-{dataset}-{shard}.txt")


if __name__ == "__main__":
    restore(sys.argv[1], int(sys.argv[2]))
