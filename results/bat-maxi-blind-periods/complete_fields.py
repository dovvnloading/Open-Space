"""Complete post-selection field checks using MAXI-only neighbours from its own roster.

The original confusion.jsonl is retained. No additional blind-frequency tests are
performed. Unlisted sources and imperfect sky-response deblending remain limits.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from field_checks import field_members, inspect_field, load_archives, read_rows
from ingest import fetch, load_curve


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_npz(path):
    with np.load(path, allow_pickle=False) as d:
        return tuple(d[k] for k in ("t", "y", "dy"))


def missing_ids(row):
    return {n["id"] for f in row["fields"] if f["band"] == "maxi"
            for n in f.get("neighbors", [])
            if n["status"] == "unresolved missing neighbor curve" and not n.get("curve_key")}


def missing_product(ident, maxi, out):
    obj = maxi[ident]
    record = {"id":ident, "aliases":obj["aliases"], "position":obj["position"], "page_url":obj["page"]}
    try:
        _, record["page"] = fetch(obj["page"])
    except Exception as exc:
        record["page_error"] = str(exc)[:300]
    row = {"id":ident, "maxi_id":ident, "maxi_page":obj["page"], "maxi_pages":obj["pages"]}
    record["product"] = load_curve(row, "maxi", out)
    if record["product"]["status"] == "ok":
        path = out/"curves"/f"{ident}.maxi.npz"
        record["frozen_curve"] = str(path.relative_to(out))
        record["frozen_sha256"] = sha(path)
    return record


def complete(evidence, out, snapshot, frozen):
    source = evidence/"validation"/"confusion.jsonl"
    original = list(read_rows(source))
    inv = json.loads((evidence/"inventory.json").read_text())
    bat = {b["id"]:b for b in inv["sources"]}
    maxi = {m["id"]:m for m in inv["maxi_sources"]}
    inverse = {}
    for b in bat.values():
        if b.get("maxi_id"):
            inverse.setdefault(b["maxi_id"], []).append(b["id"])
    (out/"curves").mkdir(parents=True, exist_ok=True)
    ids = sorted({s for r in original for s in missing_ids(r)})
    records = {}
    if frozen:
        stored = json.loads((frozen/"field-products.json").read_text())["sources"]
    for ident in ids:
        if frozen:
            entry = stored[ident]
            if entry.get("frozen_curve"):
                p = Path(entry["frozen_curve"])
                path = frozen/p
                if path.parent != frozen/"curves" or sha(path) != entry["frozen_sha256"]:
                    raise ValueError(f"unexpected path or checksum in frozen input {ident}")
                if path.resolve() != (out/p).resolve():
                    (out/p).write_bytes(path.read_bytes())
        else:
            entry = missing_product(ident, maxi, out)
        records[ident] = entry
    (out/"field-products.json").write_text(json.dumps({"original_field_sha256":sha(source),
        "purpose":"Post-selection confusion check at the daily survivor frequency, using the original field classifier.",
        "sources":records}, indent=2, allow_nan=False)+"\n")
    entries = []
    wanted = set()
    for row in original:
        if not missing_ids(row):
            entries.append(None)
            continue
        fields = field_members(row, bat, maxi, inverse)
        for f in fields:
            if f.get("target_curve"):
                wanted.add(f["target_curve"])
            for n in f["neighbors"]:
                if n.get("curve_key"):
                    wanted.add(n["curve_key"])
                elif f["band"]=="maxi" and n["id"] in records:
                    n["curve_key"] = "maxi-only/"+n["id"]
        entries.append(fields)
    curves = load_archives(sorted(snapshot.glob("curves-*.tar.gz")), wanted)
    for ident, r in records.items():
        if r.get("frozen_curve"):
            curves["maxi-only/"+ident] = load_npz(out/r["frozen_curve"])
    with (out/"confusion-completed.jsonl").open("w") as handle:
        for row, fields in zip(original, entries):
            result = inspect_field(row, fields, curves) if fields is not None else row
            handle.write(json.dumps(result, separators=(",", ":"), allow_nan=False)+"\n")
            if fields is not None:
                for f in result["fields"]:
                    for n in f.get("neighbors", []):
                        if f["band"]=="maxi" and n["id"] in records:
                            print(row["source_id"], row["period_days"], n["id"], n["status"],
                                  "amplitude ratio", n.get("amplitude_ratio_neighbor_to_target"),
                                  "phase offset", n.get("phase_offset_cycles"))
    print("Completed", len(original), "post-selection fields;", len(ids), "additional MAXI-only products")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--snapshot", type=Path, help="Folder containing the original curve shards")
    p.add_argument("--frozen", type=Path, help="Replay additional neighbours from a previously written output folder")
    a = p.parse_args()
    complete(a.evidence, a.output, a.snapshot or a.evidence, a.frozen)
