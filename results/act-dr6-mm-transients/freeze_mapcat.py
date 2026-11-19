#!/usr/bin/env python3
"""Freeze public ACT depth-1 MapCat SQLite into immutable rows with SHA-256 receipts.

INVENTORY ONLY. Does not download FITS data, find transients or constrain rates.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
import urllib.request

ROOT = "https://object-arbutus.alliancecan.ca/f620008d8888477e9fc9e5dca514fbc3:sodacan-public/maps/act"
COLUMNS = ("map_id", "map_name", "map_path", "ivar_path", "mean_time_path", "frequency", "tube_slot", "ctime", "start_time", "stop_time")
MAX_CATALOG_BYTES = 300 * 1024 * 1024


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_catalog(target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(ROOT + "/mapcat.db", headers={"User-Agent": "Open-Space-ACT-inventory/1.0"})
    with urllib.request.urlopen(request, timeout=90) as reply, target.open("wb") as out:
        content_length = int(reply.headers.get("Content-Length", "0"))
        if content_length > MAX_CATALOG_BYTES:
            raise ValueError("MapCat exceeds 300 MiB safety limit")
        total = 0
        while data := reply.read(1024 * 1024):
            total += len(data)
            if total > MAX_CATALOG_BYTES:
                raise ValueError("MapCat exceeded 300 MiB safety limit during transfer")
            out.write(data)
    with target.open("rb") as stream:
        header = stream.read(16)
    if target.stat().st_size < 100 or header != b"SQLite format 3\x00":
        raise ValueError("Downloaded MapCat is not a valid SQLite header")


def freeze(db: Path, out: Path) -> dict:
    out.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(f"file:{db.resolve()}?mode=ro", uri=True) as conn:
        conn.execute("PRAGMA query_only=ON")
        present = {r[1] for r in conn.execute('PRAGMA table_info("depth_one_maps")')}
        missing = set(COLUMNS) - present
        if missing:
            raise ValueError(f"Unsupported MapCat schema; missing columns: {sorted(missing)}")
        sql = "SELECT " + ", ".join('"' + c + '"' for c in COLUMNS) + ' FROM "depth_one_maps" ORDER BY "map_name", "map_id"'
        tally = Counter()
        empty = Counter()
        n = 0
        last_name = None
        with out.open("w", encoding="utf-8", newline="\n") as handle:
            for row in conn.execute(sql):
                obj = dict(zip(COLUMNS, row))
                if last_name is not None and obj["map_name"] == last_name:
                    raise ValueError("Duplicate map_name: " + str(last_name))
                last_name = obj["map_name"]
                for field in ("map_path", "ivar_path"):
                    if not obj[field]:
                        empty[field] += 1
                tally[str(obj["frequency"])] += 1
                handle.write(json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str) + "\n")
                n += 1
    if n == 0:
        raise ValueError("No depth_one_maps: no survey inventory can be frozen")
    return {
        "status": "catalog_inventory_only_not_FITS_verified",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": ROOT + "/mapcat.db",
        "source_sha256": digest(db),
        "source_bytes": db.stat().st_size,
        "manifest_sha256": digest(out),
        "manifest_rows": n,
        "frequency_counts_as_catalogued": dict(sorted(tally.items())),
        "missing_required_product_paths": dict(empty),
        "sql": sql,
        "actual_fits_downloads": 0,
        "actual_fits_sha256_computed": 0,
        "observational_search_completed": False,
        "rate_limit_valid": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=Path("mapcat.db"))
    parser.add_argument("--out", type=Path, default=Path("act-depth1-manifest.jsonl"))
    parser.add_argument("--receipt", type=Path, default=Path("act-depth1-inventory-receipt.json"))
    parser.add_argument("--download", action="store_true", help="Download the public SQLite catalog (not FITS maps)")
    args = parser.parse_args()
    if args.download:
        download_catalog(args.db)
    if not args.db.is_file():
        parser.error("No MapCat database at --db (pass --download or stage the SQLite file)")
    rec = freeze(args.db, args.out)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(rec, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
