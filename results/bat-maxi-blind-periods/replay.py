"""Re-run a frozen science shard, without contacting an archive. Python >=3.12."""
from __future__ import annotations

import argparse
import gzip
import json
import math
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path


def equal_values(actual,expected):
    if type(actual) is not type(expected):
        return False
    if isinstance(actual,dict):
        return actual.keys()==expected.keys() and all(equal_values(actual[k],expected[k]) for k in actual)
    if isinstance(actual,list):
        return len(actual)==len(expected) and all(equal_values(x,y) for x,y in zip(actual,expected))
    if isinstance(actual,float):
        return math.isclose(actual,expected,rel_tol=1e-8,abs_tol=2e-12)
    return actual==expected


def verify(shard,baseline,reproduced,receipt):
    name=f"search-{shard:02}.jsonl.gz"
    with gzip.open(baseline/name,"rt") as f:
        before=[json.loads(line) for line in f]
    with gzip.open(reproduced/name,"rt") as f:
        after=[json.loads(line) for line in f]
    if len(before)!=len(after):
        raise ValueError("Replay source count changed")
    mismatches=[a.get("roster_id") for a,b in zip(before,after) if not equal_values(a,b)]
    result=dict(shard=shard,records=len(before),mismatched_ids=mismatches,float_relative_tolerance=1e-8,
                float_absolute_tolerance=2e-12,comparison="All JSON fields after decompression")
    if receipt is not None:
        receipt.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result))
    if mismatches:
        raise ValueError("Replay differed from the frozen numerical outputs")


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shard", type=int)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--snapshot", type=Path, default=root / "evidence")
    parser.add_argument("--verify",action="store_true",help="Compare all fields with the committed numerical results.")
    parser.add_argument("--receipt",type=Path,help="Write the verification outcome as JSON.")
    args = parser.parse_args()
    archives = sorted(args.snapshot.glob("curves-*.tar.gz"))
    if not archives:
        raise SystemExit("The snapshot has no light-curve archives.")
    with tempfile.TemporaryDirectory() as directory:
        with tarfile.open(args.snapshot / f"curves-{args.shard:02}.tar.gz") as tar:
            for member in tar:
                if not member.isfile() or Path(member.name).parent != Path("curves") or not member.name.endswith(".npz"):
                    raise ValueError("Unexpected member path in curve snapshot")
                tar.extract(member, path=directory, filter="data")
        receipt = args.snapshot / f"ingest-{args.shard}.jsonl"
        if receipt.is_file():
            shutil.copy(receipt, directory)
        subprocess.run([
            sys.executable, str(root / "run_shard.py"),
            "--inventory", str(args.snapshot / "inventory.json"), "--work", directory,
            "--out", str(args.destination), "--shard", str(args.shard), "--shards", str(len(archives))
        ], check=True)
    if args.verify:
        verify(args.shard,args.snapshot,args.destination,args.receipt)


if __name__ == "__main__":
    main()
