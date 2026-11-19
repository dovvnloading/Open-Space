#!/usr/bin/env python3
"""Audit frozen depth-1 catalog metadata without opening science FITS.

The distinct map-name stem counts scan groups, NOT independent sky visits.
Duration between catalog start and stop is *not* time on an individual pixel.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime
import gzip
import json
from pathlib import Path
import statistics

from freeze_mapcat import digest

EXPECTED_SHA256 = "b6a0c020e3b03bc23aac66536ae1b31cb8889cbc3295c1d2acf16f8de3952503"


def audit(path: Path) -> dict:
    if digest(path) != EXPECTED_SHA256:
        raise ValueError("Manifest changed since scope freeze")
    per_year = Counter()
    per_combo = Counter()
    group_members = defaultdict(list)
    duration_minutes = []
    n = 0
    min_time, max_time = None, None
    with gzip.open(path, "rt", encoding="utf-8") as fd:
        for line in fd:
            r = json.loads(line)
            n += 1
            name = r["map_name"]
            if not name or "_pa" not in name:
                raise ValueError(f"Unexpected map name: {name}")
            stem = name.rsplit("_pa", 1)[0]
            if not stem.endswith(("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")):
                raise ValueError(f"Unexpected scan stem: {stem}")
            start = datetime.fromisoformat(r["start_time"])
            stop = datetime.fromisoformat(r["stop_time"])
            if not start < stop:
                raise ValueError(f"Nonpositive scan interval: {name}")
            minutes = (stop - start).total_seconds() / 60
            duration_minutes.append(minutes)
            group_members[stem].append((start, stop, r["frequency"]))
            per_year[str(start.year)] += 1
            per_combo[f"{r['tube_slot']}/{r['frequency']}"] += 1
            min_time = min(min_time, start) if min_time else start
            max_time = max(max_time, stop) if max_time else stop
    if n != 27512 or len(group_members) != 5579:
        raise ValueError(f"Unexpected catalog or stem group count: {n}, {len(group_members)}")
    multiplicity = Counter(map(len, group_members.values()))
    group_duration_minutes = []
    for members in group_members.values():
        start = min(m[0] for m in members)
        stop = max(m[1] for m in members)
        group_duration_minutes.append((stop - start).total_seconds() / 60)
    return {
        "status": "catalog_metadata_only_no_FITS_or_transient_measurement",
        "eligible_gzip_sha256": EXPECTED_SHA256,
        "eligible_map_array_band_rows": n,
        "distinct_map_name_scan_stems": len(group_members),
        "not_independent_spacetime_trials": True,
        "rows_per_scan_stem": dict(sorted(multiplicity.items())),
        "catalog_rows_by_start_year": dict(sorted(per_year.items())),
        "catalog_rows_by_array_band": dict(sorted(per_combo.items())),
        "earliest_start_time_no_offset": min_time.isoformat() if min_time else None,
        "latest_stop_time_no_offset": max_time.isoformat() if max_time else None,
        "median_row_start_stop_minutes": round(statistics.median(duration_minutes), 3),
        "sum_scan_stem_time_hours_uncorrected_for_sky_overlap": round(sum(group_duration_minutes) / 60, 2),
        "pixel_exposure_minutes": None,
        "usable_pixel_footprint_deg2": None,
        "astrophysical_transient_results": None,
        "rate_limit": None,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    result = audit(args.manifest)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
