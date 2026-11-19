#!/usr/bin/env python3
"""Ensure a changed image-analysis method never silently inherits old results.

Preserve every intermediate field/summary as an audit archive, then restart
from the *same* frozen S3 visit manifest. No re-selection of sources.
"""
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent
O=P/"output"
SOURCE=P/"run.py"
MARKER=O/"analysis_revision.json"
source_sha=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
old=json.loads(MARKER.read_text()) if MARKER.exists() else {}
if old.get("run_py_sha256")==source_sha:
    print("ANALYSIS_REVISION_ALREADY_CHECKED",source_sha[:12],flush=True)
    raise SystemExit(0)
previous_id=old.get("run_py_sha256","pre_revision_marker")[:12]
archive=O/"archived_previous_analysis"/previous_id
changed=[]
for name in ("north.json","south.json","fornax.json","summary.json","referee.json",
             "external_receipts.json"):
    p=O/name
    if p.exists():
        archive.mkdir(parents=True,exist_ok=True)
        target=archive/name
        if target.exists():
            raise RuntimeError("Archive collision: "+str(target))
        p.replace(target)
        changed.append(name)
MARKER.write_text(json.dumps({
    "status":"SOURCE_HASH_PINNED_RECOMPUTE_FROM_SAME_S3_VISIT_MANIFEST",
    "run_py_sha256":source_sha,
    "prior_run_py_sha256":old.get("run_py_sha256"),
    "archived_files":changed,
    "method_revision":"DR1 DPDD bit0-invalid semantics, float64 nearest-neighbour flags, prior-North v5 compact-source calibration",
    "source_data_unchanged":True
},indent=2,sort_keys=True)+"\n")
print("ANALYSIS_VERSION",source_sha[:12],"archived",changed,flush=True)
