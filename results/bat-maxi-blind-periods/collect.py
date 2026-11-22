"""Join blind numerical survivors to identifiers before any literature search."""
import argparse
import collections
import gzip
import json
from pathlib import Path


def collect(folder):
    inventory = json.loads((folder / "inventory.json").read_text())
    rows = {row["id"]: row for row in inventory["sources"]}
    tally = collections.Counter()
    seen = set()
    survivors = []
    for path in sorted(folder.glob("search-*.jsonl.gz")):
        with gzip.open(path, "rt", encoding="utf8") as handle:
            for line in handle:
                result = json.loads(line)
                ident = result["roster_id"]
                if ident in seen or ident not in rows:
                    raise ValueError(f"Duplicate or extra source {ident}")
                seen.add(ident)
                tally[result["status"]] += 1
                row = rows[ident]
                ra, dec = row.get("position") or (None, None)
                for candidate in result.get("candidates", []):
                    survivors.append({**candidate, "source_id": ident, "source": row["name"],
                                      "ra": ra, "dec": dec, "aliases": sorted(set(
                                          row["aliases"] + [x for x in (row.get("maxi_id"), row.get("maxi_name")) if x]))})
    missing = sorted(set(rows) - seen)
    with (folder / "survivors.jsonl").open("w") as handle:
        for survivor in survivors:
            handle.write(json.dumps(survivor, separators=(",", ":")) + "\n")
    (folder / "coverage.json").write_text(json.dumps(
        dict(bat_roster=len(rows), maxi_roster=len(inventory.get("maxi_sources", [])),
             associated=sum(bool(b.get("maxi_id")) for b in rows.values()), analyzed_rows=len(seen),
             absent_rows=missing, statuses=tally, numerical_survivors=len(survivors)),
        indent=2) + "\n")
    print((folder / "coverage.json").read_text())
    return missing


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    if collect(args.folder):
        raise SystemExit("Incomplete shard collection; partial rows were retained.")
