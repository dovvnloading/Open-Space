"""Download the archived optical inputs and verify receipt SHA-256 hashes.
Use --group des for the motion-search inputs only; all also includes ~1 GB for 2024/25.
"""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

root = Path(__file__).resolve().parent
cache = root / "receipts"
cache.mkdir(exist_ok=True)
ap = argparse.ArgumentParser()
ap.add_argument("--group", choices=["des", "all"], default="all")
args = ap.parse_args()
for row in json.loads((root / "decam_inputs.json").read_text()):
    if args.group == "des" and row["label"] not in [
        "DECam_605285",
        "DECam_605571",
        "mask_605285",
        "mask_605571",
        "template_605571",
    ]:
        continue
    path = cache / (row["label"] + ".response")
    if path.exists():
        with path.open("rb") as f:
            if hashlib.file_digest(f, "sha256").hexdigest() == row["sha256"]:
                continue
    snapshot = root / "query_snapshots" / (row["label"] + ".gz")
    if snapshot.exists():
        content = gzip.decompress(snapshot.read_bytes())
        if hashlib.sha256(content).hexdigest() != row["sha256"]:
            raise ValueError(f"Corrupt snapshot {snapshot.name}")
        path.write_bytes(content)
        continue
    req = Request(row["url"])
    if "payload" in row:
        req = Request(
            row["url"],
            data=json.dumps(row["payload"]).encode(),
            headers={"Content-Type": "application/json"},
        )
    partial = path.with_suffix(".part")
    digest = hashlib.sha256()
    with urlopen(req, timeout=180) as remote, partial.open("wb") as dest:
        while block := remote.read(1 << 20):
            digest.update(block)
            dest.write(block)
    if digest.hexdigest() != row["sha256"]:
        partial.unlink()
        raise ValueError(
            f"Archive bytes changed for {row['label']}, see recorded query and hash"
        )
    partial.replace(path)
