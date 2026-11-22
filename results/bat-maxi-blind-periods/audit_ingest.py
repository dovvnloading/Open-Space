"""Summarize completeness, observing spans and rejection causes of the frozen inputs."""
import argparse
from collections import Counter
from datetime import date, timedelta
import json
from pathlib import Path
from statistics import median


def summarize(data):
    if not data:
        return None
    return {"min": min(data), "median": median(data), "max": max(data)}


def iso(mjd):
    return str(date(1858, 11, 17) + timedelta(days=int(mjd)))


def audit(folder):
    sources = json.loads((folder / "inventory.json").read_text())["sources"]
    lookup = {r["id"]: r for r in sources}
    seen = {}
    for file in sorted(folder.glob("ingest-*.jsonl")):
        for line in file.read_text().splitlines():
            row = json.loads(line)
            if row["id"] not in lookup or row["id"] in seen:
                raise ValueError(f"duplicate or unknown ingestion row {row['id']}")
            seen[row["id"]] = row
    missing = sorted(lookup.keys() - seen.keys())
    report = {"bat_roster_count": len(lookup), "receipt_rows": len(seen),
              "missing_ids": missing, "match_methods": Counter(b["match"] for b in sources)}
    for label, instrument, rows in (
            ("bat_all", "bat", sources),
            ("bat_weak", "bat", [b for b in sources if "/weak/" in b["bat_page"]]),
            ("bat_main", "bat", [b for b in sources if "/weak/" not in b["bat_page"]]),
            ("maxi_matched", "maxi", [b for b in sources if b["maxi_id"]])):
        products = [(b["id"], seen[b["id"]][instrument]) for b in rows if b["id"] in seen]
        ok = [r for _, r in products if r["status"] == "ok"]
        first = summarize([r["first_mjd"] for r in ok])
        last = summarize([r["last_mjd"] for r in ok])
        report[label] = {
            "roster": len(rows), "status": Counter(r["status"] for _,r in products),
            "filtered_rows": sum(r["rows_raw"] - r["rows_used"] for r in ok),
            "retained_rows": sum(r["rows_used"] for r in ok),
            "retained_days_per_source": summarize([r["rows_used"] for r in ok]),
            "first_mjd": first, "last_mjd": last,
            "first_date": {k:iso(v) for k,v in first.items()} if first else None,
            "last_date": {k:iso(v) for k,v in last.items()} if last else None,
            "source_span_days": summarize([r["last_mjd"]-r["first_mjd"] for r in ok]),
            "failed_ids": [s for s,r in products if r["status"] != "ok"],
            "negative_rate_days": sum(r["negative_rates"] for r in ok)
        }
    (folder/"ingest-audit.json").write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
    for k in ("bat_main", "bat_weak", "maxi_matched"):
        cohort = report[k]
        print(k, cohort["status"], cohort["first_date"], cohort["last_date"])
    if missing:
        raise ValueError(f"missing {len(missing)} archive receipts; see ingest-audit.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    audit(parser.parse_args().folder)
