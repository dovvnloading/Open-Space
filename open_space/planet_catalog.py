"""Batch known-host and period checks against NASA's confirmed planet and TOI tables."""

import argparse
import csv
import io
import json
from pathlib import Path
from urllib.parse import urlencode

from .archive_index import EXOPLANETS, record_response, retrieve

TABLES = {
    "pscomppars": ("tic_id", "pl_name", "tic_id,pl_name,hostname,pl_orbper,ra,dec"),
    "toi": ("tid", "toi", "tid,toi,tfopwg_disp,pl_orbper,ra,dec"),
}


def snapshot(targets, cache):
    ids = sorted({int(t) for t in targets})
    if not ids or min(ids) < 1:
        raise ValueError("Supply positive TIC IDs")
    matches = {str(tic): [] for tic in ids}
    receipts = []
    for start in range(0, len(ids), 100):
        batch = ids[start : start + 100]
        for table, (id_key, name_key, columns) in TABLES.items():
            vals = ",".join(
                f"'TIC {tic}'" if id_key == "tic_id" else str(tic) for tic in batch
            )
            query = f"SELECT {columns} FROM {table} WHERE {id_key} IN ({vals}) ORDER BY {id_key},{name_key}"
            url = EXOPLANETS + "?" + urlencode({"query": query, "format": "csv"})
            raw, path = retrieve(url, cache, "text/csv")
            reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
            if set(reader.fieldnames or []) != set(columns.split(",")):
                raise ValueError(
                    f"Unexpected {table} response; catalog lookup incomplete"
                )
            for row in reader:
                tic = int(row[id_key].removeprefix("TIC ").strip())
                if tic not in batch:
                    raise ValueError("Unexpected target in catalog result")
                matches[str(tic)].append(
                    {
                        "table": table,
                        "object": row[name_key],
                        "archive_disposition": row.get("tfopwg_disp") or None,
                        "period_days": (
                            float(row["pl_orbper"]) if row["pl_orbper"] else None
                        ),
                        "ra_deg": row["ra"],
                        "dec_deg": row["dec"],
                    }
                )
            receipts.append(
                {**record_response(raw, path, url), "table": table, "query": query}
            )
    return {
        "schema_version": 1,
        "match_method": "exact TIC identifier",
        "queries": receipts,
        "targets": matches,
    }


def assess(tic, period, catalog):
    if str(tic) not in catalog["targets"]:
        return {"status": "not_queried", "matches": []}
    rows = []
    for target in catalog["targets"][str(tic)]:
        reference = target["period_days"]
        comparisons = (
            [
                (abs(period / (reference * factor) - 1), label)
                for factor, label in [
                    (1, "period_match"),
                    (0.5, "half_period"),
                    (2, "double_period"),
                ]
            ]
            if reference and period
            else []
        )
        delta, nearest = min(comparisons, default=(None, None))
        rows.append(
            {
                **target,
                "period_relation": (
                    nearest
                    if delta is not None and delta <= 0.02
                    else "no_close_period"
                ),
                "nearest_fractional_offset": delta,
            }
        )
    associated = any(row["period_relation"] != "no_close_period" for row in rows)
    status = (
        "catalog_ephemeris_match"
        if associated
        else (
            "catalog_host_other_or_missing_period"
            if rows
            else "no_exact_tic_entry_in_snapshot"
        )
    )
    return {"status": status, "matches": rows}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("plan", type=Path)
    p.add_argument("--cache", type=Path, default=Path("runs/catalog-cache"))
    args = p.parse_args()
    plan = json.loads(args.plan.read_bytes())
    result = snapshot([r["tic_id"] for r in plan["targets"]], args.cache)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
