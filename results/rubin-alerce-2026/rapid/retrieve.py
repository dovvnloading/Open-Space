"""Replay the frozen selection, light-curve and nearby-source queries."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    args.destination.mkdir(parents=True, exist_ok=True)
    for item in json.loads((ROOT / "data/queries.json").read_text()):
        params = dict(
            REQUEST="doQuery",
            LANG="ADQL",
            FORMAT="csv",
            MAXREC=100000,
            QUERY=item["adql"],
        )
        with urlopen(
            "https://tap.alerce.online/tap/sync?" + urlencode(params), timeout=180
        ) as response:
            body = response.read()
        records = csv.DictReader(io.StringIO(body.decode()))
        if "oid" not in (records.fieldnames or []):
            raise ValueError(f"Unexpected response for {item['key']}")
        count = sum(1 for row in records)
        digest = hashlib.sha256(body).hexdigest()
        (args.destination / (item["key"] + ".csv")).write_bytes(body)
        print(
            item["key"],
            count,
            "same bytes" if digest == item["sha256"] else "changed bytes",
            digest,
        )
