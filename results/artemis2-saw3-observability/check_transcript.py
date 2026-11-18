"""Compare the two PCD transcript time columns against the product label."""

from argparse import ArgumentParser
import csv
from datetime import datetime, timedelta
from hashlib import sha256
import json
from pathlib import Path
from statistics import median
from xml.etree import ElementTree


CSV_SHA = "6d555f49ea7f6ea717c7d15a2cf41ef16e58cce5d55f73e781334d145d6acfdf"
LABEL_SHA = "8a88d73640366df2725111dd1730fbfbf942b2e31869160089dc47ed66680cd5"
NS = {"pds": "http://pds.nasa.gov/pds4/pds/v1"}
BASE = (
    "https://pds-imaging.jpl.nasa.gov/archive/archive/artemis/artemis2/"
    "artemis2_mission_audio/data_transcripts_pcd/fd06/"
)


def relative_offset(value):
    hours, minutes, seconds = map(int, value.split(":"))
    return timedelta(hours=hours, minutes=minutes, seconds=seconds)


def audit(table, label):
    if sha256(table.read_bytes()).hexdigest() != CSV_SHA or sha256(label.read_bytes()).hexdigest() != LABEL_SHA:
        raise ValueError("Input differs from the PDS transcript or label used in this audit")
    xml = ElementTree.parse(label)
    start = datetime.fromisoformat(xml.findtext(
        ".//pds:Time_Coordinates/pds:start_date_time", namespaces=NS).replace("Z", "+00:00"))
    offsets, candidates = [], []
    with table.open(encoding="utf-8-sig", newline="") as stream:
        for number, row in enumerate(csv.DictReader(stream), start=1):
            relative = start + relative_offset(row["audio_time"])
            assigned = datetime.fromisoformat(row["utc"]).replace(tzinfo=start.tzinfo)
            offset = (assigned - relative).total_seconds()
            offsets.append(offset)
            if number in (159, 272):
                candidates.append({
                    "record_number": number, "relative_audio_time": row["audio_time"],
                    "assigned_utc": assigned.isoformat(), "utc_from_relative": relative.isoformat(),
                    "assigned_minus_relative_seconds": offset,
                })
    return {
        "source": {"file": table.name, "raw_sha256": CSV_SHA, "url": BASE + table.name},
        "label": {"file": label.name, "raw_sha256": LABEL_SHA, "url": BASE + label.name},
        "start_utc": start.isoformat(), "transcript_rows": len(offsets),
        "offset_limits_seconds": [min(offsets), max(offsets)], "offset_median_seconds": median(offsets),
        "outside_one_second": sum(abs(offset) > 1 for offset in offsets), "cue_rows": candidates,
    }


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("pcd_csv", type=Path)
    parser.add_argument("pcd_label", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = audit(args.pcd_csv, args.pcd_label)
    with args.output.open("x", encoding="utf-8", newline="\n") as out:
        out.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
