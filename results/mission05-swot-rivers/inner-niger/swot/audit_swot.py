"""Audit reach coverage, source lineage, sampling, and SoS estimator availability.

Run: python audit_swot.py DIRECTORY_WITH_MANIFEST OUTPUT_DIRECTORY
Results are descriptive and summarize distinct reach-pass rows. Adjacent reaches
often share the same overpass; these counts are not independent sample sizes.
"""

import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from datetime import datetime
from math import isfinite
from pathlib import Path
from statistics import median

# DAWG lists neoBAM as the sixth FLPE member, currently not exposed by Hydrocron.
# LakeFlow is another exposed estimate for connected lakes, audited separately.
FLPE_NAMES = ["hivdi", "metroman", "momma", "sad", "sic4dvar"]
FLPE_FIELDS = [f"sos_{name}_q" for name in FLPE_NAMES]
EXPOSED_FIELDS = FLPE_FIELDS + ["sos_lakeflow_q"]


def number(row, field):
    try:
        value = float(row[field])
    except (TypeError, ValueError):
        return None
    return value if isfinite(value) and value != -999999999999.0 else None


def principal(row):
    return (
        int(row["reach_q"]) <= 1
        and row["partial_f"] == "0"
        and float(row["obs_frac_n"]) >= 0.8
        and int(row["xovr_cal_q"]) in (0, 1)
    )


def stricter(row):
    # Exclude suspect geolocation, the near-nadir gap, and nonzero crossover flags.
    # Distances are additional diagnostic sensitivity choices in metres.
    return (
        principal(row)
        and int(row["reach_q_b"]) & 4 == 0
        and row["xovr_cal_q"] == "0"
        and float(row["obs_frac_n"]) >= 0.95
        and number(row, "loc_offset") is not None
        and abs(float(row["loc_offset"])) <= 100
        and number(row, "xtrk_dist") is not None
        and 10000 <= abs(float(row["xtrk_dist"])) <= 60000
    )


def summarized(rows, reach, sample):
    out = {"reach_id": reach, "sample": sample, "n_rows": len(rows)}
    for q in range(4):
        out[f"n_reach_q_{q}"] = sum(int(row["reach_q"]) == q for row in rows)
    for field in ["wse", "slope", "width", "sos_consensus_q"] + EXPOSED_FIELDS:
        out[f"n_{field}"] = sum(number(row, field) is not None for row in rows)
    has_q = [row for row in rows if number(row, "sos_consensus_q") is not None]
    for name, field in zip(FLPE_NAMES, FLPE_FIELDS):
        out[f"n_consensus_equals_{name}"] = sum(
            number(row, "sos_consensus_q") == number(row, field) for row in has_q
        )
    support = Counter(sum(number(row, f) is not None for f in FLPE_FIELDS) for row in has_q)
    for n in range(len(FLPE_NAMES) + 1):
        out[f"n_consensus_with_{n}_exposed_flpe_present"] = support[n]
    for field in ["wse_u", "slope_u", "width_u"]:
        valid = [number(row, field) for row in rows if number(row, field) is not None]
        out[f"median_reported_{field}"] = median(valid) if valid else ""
    valid_times = sorted(datetime.fromisoformat(row["time_str"]) for row in rows if row["time_str"] != "no_data")
    out["first_measurement_utc"] = valid_times[0].isoformat() if valid_times else ""
    out["last_measurement_utc"] = valid_times[-1].isoformat() if valid_times else ""
    gaps = [(b-a).total_seconds()/86400 for a,b in zip(valid_times, valid_times[1:])]
    out["longest_gap_days"] = round(max(gaps), 2) if gaps else ""
    out["n_cycle_pass"] = len({(row["cycle_id"], row["pass_id"]) for row in rows})
    out["n_cycles"] = len({row["cycle_id"] for row in rows})
    for crid in ("PGC0", "PIC0", "PIC2"):
        out[f"n_{crid}"] = sum(row["crid"] == crid for row in rows)
    return out


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.inputs / "manifest.json").read_text())
    results = []
    for item in manifest["reach_requests"]:
        for key in ["body", "decoded_csv"]:
            blob = (args.inputs / item[key]["path"]).read_bytes()
            if hashlib.sha256(blob).hexdigest() != item[key]["sha256"]:
                raise ValueError(f"Hash mismatch on {item[key]['path']}")
        body = json.loads((args.inputs / item["body"]["path"]).read_bytes())
        csv_bytes = body["results"]["csv"].encode("utf-8")
        if csv_bytes != (args.inputs / item["decoded_csv"]["path"]).read_bytes():
            raise ValueError("Decoded CSV does not correspond to raw JSON")
        records = list(csv.DictReader(io.StringIO(csv_bytes.decode("utf-8"))))
        reach = item["reach"]
        for sample, rows in [
            ("all", records),
            ("primary", [r for r in records if principal(r)]),
            ("strict_geolocation", [r for r in records if stricter(r)]),
        ]:
            results.append(summarized(rows, reach, sample))
    write_csv(args.output / "swot_coverage.csv", results)


if __name__ == "__main__":
    main()
