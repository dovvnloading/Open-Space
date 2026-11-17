"""Restore the original IRSA responses from the compressed snapshots and check their identity."""
import gzip
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
for record in json.loads((ROOT/"frozen_inputs.json").read_text()):
    payload=(ROOT/record["file"]).read_bytes()
    if hashlib.sha256(payload).hexdigest()!=record["gzip_sha256"]:raise ValueError(record["file"])
    raw=gzip.decompress(payload)
    if hashlib.sha256(raw).hexdigest()!=record["sha256"]:raise ValueError(record["original"])
    (ROOT/"inputs"/record["original"]).write_bytes(raw)
