"""Check saved evidence offline, including timeline order, dates, and receipt digests."""
import csv
import hashlib
import json
import re
from datetime import datetime, timedelta
from pathlib import Path

OUT = Path(__file__).resolve().parent / "outputs"


def main():
    inventory = json.loads((OUT / "snapshot.json").read_text())
    for name, digest in inventory["products_sha256"].items():
        assert hashlib.sha256((OUT / name).read_bytes()).hexdigest() == digest, name
    requests = [json.loads(line) for line in (OUT / "inputs.jsonl").read_text().splitlines()]
    assert all(re.fullmatch("[a-f0-9]{64}", r["sha256"]) for r in requests if "sha256" in r)
    assert all("sha256" in r or "error" in r for r in requests)
    for r in requests:
        if r["method"] == "git input":
            assert hashlib.sha256(Path(r["path"]).read_bytes()).hexdigest() == r["sha256"], r["path"]
    with (OUT / "timeline.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert rows and len(rows) == inventory["timeline_rows"]
    times = [float(r["mjd"]) for r in rows if r["mjd"]]
    assert times == sorted(times)
    for row in rows:
        if row["mjd"]:
            date = datetime(1858,11,17) + timedelta(days=float(row["mjd"]))
            assert abs((date - datetime.fromisoformat(row["utc"])).total_seconds()) < 0.51
        if row["end_mjd"] and row["mjd"]:
            assert float(row["end_mjd"]) >= float(row["mjd"]), row
    print(f"Verified {len(rows)} timeline records and {len(requests)} input receipts.")


if __name__ == "__main__":
    main()
