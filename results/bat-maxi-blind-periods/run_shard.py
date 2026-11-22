"""Run the same blind search on each BAT roster entry assigned to a shard."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import traceback

import numpy as np

from search import analyze


def read_curve(root, name, instrument):
    path = root / "curves" / f"{name}.{instrument}.npz"
    if path.is_file():
        with np.load(path, allow_pickle=False) as record:
            return tuple(record[key] for key in ("t", "y", "dy"))
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--shards", type=int, required=True)
    args = parser.parse_args()
    rows = json.loads(args.inventory.read_text())["sources"]
    subset = rows[args.shard::args.shards]
    args.out.mkdir(exist_ok=True, parents=True)
    with gzip.open(args.out / f"search-{args.shard:02}.jsonl.gz", "wt", encoding="utf8") as stream:
        for index, row in enumerate(subset):
            source = row["id"]
            try:
                result = analyze(source, read_curve(args.work, source, "bat"),
                                 read_curve(args.work, source, "maxi"), n_sources=len(rows))
            except Exception:
                result = {"status": "analysis_error", "traceback": traceback.format_exc()}
            result["roster_id"] = source
            stream.write(json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n")
            stream.flush()
            if index % 10 == 0:
                print(f"Shard {args.shard}: {index+1}/{len(subset)}; {source}", flush=True)
    print(f"Completed {len(subset)} roster entries.", flush=True)


if __name__ == "__main__":
    main()
