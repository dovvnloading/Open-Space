"""Recompute statistics offline using the immutable normalized light curves."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tarfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("snapshot", type=Path)
parser.add_argument("destination", type=Path)
args = parser.parse_args()
args.destination.mkdir(parents=True, exist_ok=False)
study = Path(__file__).resolve().parent
index = args.snapshot / "index/catalog.json"
roster = json.loads(index.read_text())
owner = {row["name"]: i % 24 for i, row in enumerate(roster)}
receipts = [[] for _ in range(24)]
for line in (args.snapshot / "inputs.jsonl").read_text().splitlines():
    receipts[owner[json.loads(line)["name"]]].append(line)
for i in range(24):
    data, part = args.destination / f"input-{i}", args.destination / f"part-{i}"
    with tarfile.open(args.snapshot / f"snapshot/curves-{i}.tar.gz") as tar:
        tar.extractall(data, filter="data")
    subprocess.run([sys.executable, str(study/"survey.py"), "shard", "--index", str(index),
                    "--shard", str(i), "--shards", "24", "--data", str(data), "--out", str(part)], check=True)
    (part / "manifest.jsonl").write_text("\n".join(receipts[i]) + "\n")
subprocess.run([sys.executable, str(study/"survey.py"), "reduce", "--index", str(index),
                "--work", str(args.destination), "--out", str(args.destination/"outputs")], check=True)
