"""Seal the per-source evidence without changing the scientific disposition."""
import csv
import hashlib
import json
import os
from pathlib import Path

OUT = Path(__file__).resolve().parent / "outputs"
with (OUT / "timeline.csv").open() as f:
    rows = sum(1 for _ in csv.DictReader(f))
files = sorted(p for p in OUT.rglob("*") if p.is_file() and p.name != "snapshot.json")
(OUT / "snapshot.json").write_text(json.dumps({
    "collection_revision": os.environ.get("GITHUB_SHA"),
    "workflow_run": os.environ.get("GITHUB_RUN_ID"), "timeline_rows": rows,
    "products_sha256": {p.relative_to(OUT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
}, indent=2)+"\n")
