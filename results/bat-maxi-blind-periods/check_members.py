"""Verify all sixteen frozen numerical archives against their source-specific receipts."""
from collections import Counter
import json
from pathlib import Path
import sys
import tarfile


def verify(folder, path):
    inventory = json.loads((folder / "inventory.json").read_text())["sources"]
    archives = sorted(folder.glob("curves-*.tar.gz"))
    if not archives:
        raise ValueError("No frozen archives")
    report = []
    for i,archive in enumerate(archives):
        receipt = folder/f"ingest-{i}.jsonl"
        rows = {r["id"]:r for r in map(json.loads,receipt.read_text().splitlines())}
        roster_ids = {r["id"] for r in inventory[i::len(archives)]}
        expected = {f"curves/{source}.{band}.npz" for source,r in rows.items() for band in ("bat","maxi") if r[band]["status"]=="ok"}
        with tarfile.open(archive,"r:gz") as tf:
            names = [m.name for m in tf.getmembers() if m.isfile()]
        got = set(names)
        report.append(dict(shard=i,expected=len(expected),files=len(names),missing_files=sorted(expected-got),
                           extra_files=sorted(got-expected),duplicate_members=sorted(k for k,v in Counter(names).items() if v>1),
                           missing_source_rows=sorted(roster_ids-rows.keys()),extra_source_rows=sorted(rows.keys()-roster_ids)))
    path.write_text(json.dumps(report,indent=2)+"\n")
    if any(r[k] for r in report for k in ("missing_files","extra_files","duplicate_members","missing_source_rows","extra_source_rows")):
        raise RuntimeError("The archive member checks did not reproduce the complete frozen roster.")
    print("Matched every expected light-curve member to its shard and instrument.")


if __name__=="__main__":
    verify(Path(sys.argv[1]),Path(sys.argv[2]))
