#!/usr/bin/env python3
"""Materialize the six ACT DR6 depth-1 map/ivar inventory, excluding pa7 f030/f040.

Filters ONLY catalog metadata; inclusion is not verified image/ivar coverage.
Compresses a fixed-order, timestamp-free JSONL with an immutable checksum.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import json
from pathlib import Path

from freeze_mapcat import digest

EXPECTED_SOURCE_SHA256 = "c4e73ee1a0c455dddbc7c01160a08b37cd5264020096792d1ba2beedd53d8efd"
EXPECTED_MANIFEST_SHA256 = "a78aa7ebdb263153f1a3a86a5ad2a67f93b924c11ebe0267e54a7e3f16df8f85"
COMBOS = {
    ("pa5", "f090"), ("pa6", "f090"),
    ("pa4", "f150"), ("pa5", "f150"), ("pa6", "f150"),
    ("pa4", "f220"),
}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--catalog-receipt", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--receipt", type=Path, required=True)
    a = p.parse_args()
    receipt = json.loads(a.catalog_receipt.read_text())
    if receipt.get("source_sha256") != EXPECTED_SOURCE_SHA256:
        p.error("Public MapCat bytes changed; this is a NEW dataset, do not silently replace freeze")
    if digest(a.manifest) != EXPECTED_MANIFEST_SHA256 or receipt.get("manifest_sha256") != EXPECTED_MANIFEST_SHA256:
        p.error("Full catalog manifest differs from frozen inventory")
    counts = Counter()
    rejected = Counter()
    nall = 0
    nsel = 0
    last_name = ""
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.manifest.open("r", encoding="utf-8") as source, a.out.open("wb") as handle:
        with gzip.GzipFile(fileobj=handle, mode="wb", filename="", mtime=0) as target:
            for line in source:
                nall += 1
                row = json.loads(line)
                name = row["map_name"]
                if name <= last_name:
                    raise ValueError("MapCat name order changed or duplicate name")
                last_name = name
                combo = (row["tube_slot"], row["frequency"])
                if combo not in COMBOS or not ("2017-" <= row["start_time"] < "2023-"):
                    rejected[str(combo)] += 1
                    continue
                if not row.get("map_path") or not row.get("ivar_path") or not row.get("mean_time_path"):
                    raise ValueError("Eligible map missing map/ivar/time path")
                if "_" + combo[0] + "_" + combo[1] not in name:
                    raise ValueError("Catalog map_name disagrees with tube/frequency")
                target.write((json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode())
                counts[str(combo)] += 1
                nsel += 1
    if nall != 31261 or nsel != 27512:
        raise ValueError(f"Unexpected filtered denominator {nall} -> {nsel}; stop")
    summary = {
        "type": "frozen_catalog_eligibility_only",
        "scientific_result": False,
        "full_catalog_rows": nall,
        "eligible_depth_one_map_array_rows": nsel,
        "excluded_rows": nall - nsel,
        "excluded_by_band_array": dict(sorted(rejected.items())),
        "eligible_by_band_array": dict(sorted(counts.items())),
        "source_mapcat_sha256": EXPECTED_SOURCE_SHA256,
        "full_manifest_sha256": EXPECTED_MANIFEST_SHA256,
        "eligible_gzip_sha256": digest(a.out),
        "eligible_gzip_bytes": a.out.stat().st_size,
        "row_fields": ["map_id", "map_name", "map_path", "ivar_path", "mean_time_path",
                       "frequency", "tube_slot", "ctime", "start_time", "stop_time"],
        "fits_downloads": 0,
        "matched_filter_runs": 0,
        "independent_sky_visits_derived": False,
        "effective_survey_area_time": None,
        "rate_limit": None,
    }
    a.receipt.parent.mkdir(parents=True, exist_ok=True)
    a.receipt.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
