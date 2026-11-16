"""Check the topology and coincident principal-screen Q at the Mopti confluence.

Usage: python closure.py DIRECTORY_WITH_NINE_REACH_CSVS
"""
import csv
import io
import json
import statistics
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
WEST, EAST, OUTLET = ("14291000011", "14281000011", "14279000471")


def linked(row, field):
    return {v.strip() for v in row[field].split(",") if v.strip() != "no_data"}


def topology():
    for reach in (WEST, EAST, OUTLET):
        payload = json.loads((ROOT/f"{reach}_topology.json").read_bytes())
        for r in csv.DictReader(io.StringIO(payload["results"]["csv"])):
            if reach == OUTLET:
                assert linked(r, "rch_id_up") == {WEST, EAST}
            else:
                assert linked(r, "rch_id_dn") == {OUTLET}


def series(directory, reach):
    with (directory/f"{reach}.csv").open() as f:
        records = [
            row for row in csv.DictReader(f)
            if int(row["reach_q"]) <= 1 and row["partial_f"] == "0"
            and float(row["obs_frac_n"]) >= 0.80 and int(row["xovr_cal_q"]) <= 1
        ]
    return {(r["cycle_id"], r["pass_id"]): r for r in records}


def main():
    topology()
    directory = Path(sys.argv[1])
    western, eastern, downstream = [series(directory, r) for r in [WEST, EAST, OUTLET]]
    joined = sorted(western.keys() & eastern.keys() & downstream.keys())
    result = []
    for cycle, pass_id in joined:
        key = (cycle, pass_id)
        rows = [group[key] for group in (western, eastern, downstream)]
        q = [float(r["sos_consensus_q"]) for r in rows]
        if min(q) <= 0:
            continue
        times = [datetime.fromisoformat(r["time_str"]) for r in rows]
        difference = q[0] + q[1] - q[2]
        result.append({
            "cycle_id": cycle, "pass_id": pass_id, "downstream_utc": rows[2]["time_str"],
            "sample_spread_seconds": (max(times)-min(times)).total_seconds(),
            "west_q_m3s": q[0], "east_q_m3s": q[1], "downstream_q_m3s": q[2],
            "arms_minus_downstream_m3s": difference,
            "fraction_of_downstream": difference/q[2],
            "common_granule": rows[2]["granuleUR"],
        })
        assert len({r["granuleUR"] for r in rows}) == 1
    with (ROOT/"confluence_closure.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, result[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(result)
    frac = [r["fraction_of_downstream"] for r in result]
    summary = {
        "complete_principal_qc_orbits_without_discharge_condition":len(joined),
        "positive_L4_consensus_triples":len(result),
        "max_timestamp_separation_s":max(r["sample_spread_seconds"] for r in result),
        "median_signed_fraction":statistics.median(frac),
        "median_absolute_fraction":statistics.median(map(abs,frac)),
        "range_fraction":[min(frac),max(frac)],
        "by_orbit":dict(Counter(r["pass_id"] for r in result)),
        "cycles":len({r["cycle_id"] for r in result}),
        "first_date":result[0]["downstream_utc"],"last_date":result[-1]["downstream_utc"],
    }
    (ROOT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__ == "__main__":
    main()
