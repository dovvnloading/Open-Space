"""Re-run a frozen science shard, without contacting an archive. Python >=3.12."""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shard", type=int)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--snapshot", type=Path, default=root / "evidence")
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


if __name__ == "__main__":
    main()
