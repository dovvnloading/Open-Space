"""Freeze brightness pairs with recorded ASPIICS v03 continuum lineage."""

from argparse import ArgumentParser
from collections import Counter
from datetime import datetime
import gzip
import hashlib
import json
from pathlib import Path


BASE = "https://p3sc.oma.be/datarepfiles/"
GROUPS = ((795, 2), (772, 0))
SAMPLE_SIZE = 5
MAX_OFFSET = 90.0


def epoch(row):
    return datetime.fromisoformat(row["DATE-AVG"])


def parents(row):
    return {row[key] for key in ("PARENT1", "PARENT2", "PARENT3") if row[key] is not None}


def product(row):
    return {
        "name": row["name"], "url": BASE + row["datalocation"],
        "cycle": row["CYCLE_ID"], "date_avg": row["DATE-AVG"],
        "quality": row["QUALITY"], "n_parent": row["N_PARENT"],
        "listed_parents": sorted(parents(row)),
    }


def select(rows, quality):
    brightness = [row for row in rows if row["PROD_ID"] == "Total brightness" and row["QUALITY"] == 0]
    lines = [row for row in rows if row["PROD_ID"] == "Green line"]
    result = []
    for line in lines:
        wb = {parent for parent in parents(line) if parent.startswith("aspiics_wb_l2_")}
        donors = [row for row in brightness if wb and wb <= parents(row)]
        if len(donors) != 1:
            raise ValueError(f"{len(donors)} matching wideband parent sets for {line['name']}")

        donor = donors[0]
        offset = (epoch(donor) - epoch(line)).total_seconds()
        following = [row for row in brightness if epoch(row) > epoch(line)]
        next_brightness = min(following, key=epoch, default=None)
        next_offset = (epoch(next_brightness) - epoch(line)).total_seconds() if following else None
        result.append({
            "line": product(line),
            "donor": product(donor),
            "donor_offset_s": offset,
            "next_brightness": product(next_brightness) if following else None,
            "next_offset_s": next_offset,
            "cross_iteration": line["CYCLE_ID"] != donor["CYCLE_ID"],
        })
    eligible = [row for row in result if row["line"]["quality"] == quality
                and -MAX_OFFSET <= row["donor_offset_s"] < 0
                and row["next_offset_s"] is not None and 0 < row["next_offset_s"] <= MAX_OFFSET]
    if len(eligible) < SAMPLE_SIZE:
        raise ValueError("Too few eligible brightness brackets")
    indices = [i * (len(eligible) - 1) // (SAMPLE_SIZE - 1) for i in range(SAMPLE_SIZE)]
    return result, [eligible[i] for i in indices], len(eligible), indices


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("metadata_dir", type=Path, help="Directory with the two gzipped API receipts")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = {}
    for orbit, quality in GROUPS:
        path = args.metadata_dir / f"orbit{orbit}.json.gz"
        raw = gzip.decompress(path.read_bytes())
        rows = json.loads(raw)
        if any(row["ORBIT_ID"] != orbit for row in rows):
            raise ValueError(f"Wrong orbit in {path.name}")
        lineage, sample, eligible, indices = select(rows, quality)
        result[str(orbit)] = {
            "api_rows": len(rows), "api_sha256": hashlib.sha256(raw).hexdigest(),
            "fe_quality": dict(sorted(Counter(row["line"]["quality"] for row in lineage).items())),
            "sample_quality": quality, "eligible_brackets": eligible, "sample_indices": indices,
            "sample": sample, "lineage": lineage,
        }
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
