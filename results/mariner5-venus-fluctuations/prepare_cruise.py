"""Select quality-filtered Mariner 5 cruise frames from verified source bytes."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from statistics import fmean, pstdev


def selected(instant, rate):
    day = instant.date().isoformat()
    if ("1967-10-01" <= day <= "1967-10-15" or "1967-10-22" <= day <= "1967-11-05"):
        if rate != 0:
            return None
        seconds = instant.hour * 3600 + instant.minute * 60 + instant.second + instant.microsecond / 1e6
        return "cruise" if any(abs(seconds - hour * 3600) <= 20 * 60 for hour in (3, 9, 15, 21)) else None
    if (datetime(1967, 6, 14, 11, 55, tzinfo=timezone.utc) <= instant
            < datetime(1967, 6, 14, 12, 55, tzinfo=timezone.utc)):
        return "earth" if rate == 1 else None
    return None


def run(path, receipt, output, stats_path):
    source_receipt = json.loads(receipt.read_text(encoding="utf-8"))
    source = next(item for item in source_receipt if item["file"] == "mag_cruise_data.csv")
    expected = source["raw_sha256"]
    if not isinstance(expected, str) or len(expected) != 64:
        raise ValueError("Receipt must contain raw_sha256 or sha256")
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    actual = digest.hexdigest()
    if actual.lower() != expected.lower() or path.stat().st_size != source["bytes"]:
        raise ValueError("Cruise CSV hash does not match receipt")
    counts = {group: {"selected": 0, "rejected_tags": 0, "rejected_values": 0, "accepted": 0}
              for group in ("cruise", "earth")}
    columns = ["utc", "rate", "record_n", "sample", "mean_magnitude_nt", "vector_std_nt",
               "magnitude_std_nt", "ntype", "ij1", "ij2", "ij3", "ir1", "ir2", "ir3", "mf6", "mf10", "mf11", "mf12"]
    with path.open("r", encoding="ascii", newline="") as inp, output.open("w", encoding="utf-8", newline="") as out:
        reader = csv.DictReader(inp)
        writer = csv.DictWriter(out, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for row in reader:
            row = {key.strip(): value.strip() for key, value in row.items()}
            instant = datetime.fromisoformat(row["UTC"].replace("Z", "+00:00"))
            if instant.utcoffset() is None:
                raise ValueError("Source UTC lacks timezone")
            rate = int(row["RATE"])
            group = selected(instant, rate)
            if group is None:
                continue
            count = counts[group]
            count["selected"] += 1
            tags = {key: int(row[key.upper()]) for key in ("ntype", "ij1", "ij2", "ij3", "ir1", "ir2", "ir3", "mf6", "mf10", "mf11", "mf12")}
            if (any(tags[key] != 0 for key in ("ij1", "ij2", "ij3", "ir1", "ir2", "ir3", "mf6", "mf11", "mf12"))
                    or tags["mf10"] not in ((0, 3) if group == "earth" else (0,))):
                count["rejected_tags"] += 1
                continue
            vectors = [[float(row[f"IPM{axis}{sample}"]) for axis in (1, 2, 3)] for sample in (1, 2, 3)]
            if not all(math.isfinite(value) and abs(value) < 51.1 for vector in vectors for value in vector):
                count["rejected_values"] += 1
                continue
            magnitudes = [math.sqrt(sum(value * value for value in vector)) for vector in vectors]
            vector_std = math.sqrt(sum(pstdev([v[axis] for v in vectors]) ** 2 for axis in range(3)))
            writer.writerow({"utc": instant.isoformat().replace("+00:00", "Z"), "rate": rate,
                             "record_n": row["RECORD_N"], "sample": group,
                             "mean_magnitude_nt": format(fmean(magnitudes), ".12g"),
                             "vector_std_nt": format(vector_std, ".12g"),
                             "magnitude_std_nt": format(pstdev(magnitudes), ".12g"), **tags})
            count["accepted"] += 1
    output_sha = hashlib.sha256(output.read_bytes()).hexdigest()
    stats_path.write_text(json.dumps({"source": source["file"], "source_sha256": actual,
                                      "receipt": receipt.name, "output_sha256": output_sha,
                                      "selection": counts,
                                      "policy": "IR=IJ=MF6=MF11=MF12=0; MF10=0, or 0/3 for Earth; abs component <51.1 nT"},
                                     indent=2) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--stats", type=Path, required=True)
    args = parser.parse_args()
    run(args.input, args.receipt, args.output, args.stats)
