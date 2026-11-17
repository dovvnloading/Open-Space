"""Query ESA Q1 photometric redshifts by exact MER identifier. Stdlib only."""
import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from io import StringIO
from pathlib import Path
from time import sleep
from urllib.parse import urlencode
from urllib.request import Request, urlopen

FIELDS = ("object_id", "phz_median", "phz_70_int", "phz_95_int", "phz_mode_1",
          "phz_mode_1_area", "phz_mode_2", "phz_mode_2_area", "phz_flags", "phz_classification")
ENDPOINT = "https://eas.esac.esa.int/tap-server/tap/sync"


def ajax(ids, directory):
    query = "SELECT " + ",".join(FIELDS) + " FROM q1.phz_photo_z WHERE object_id IN (" + ",".join(ids) + ") ORDER BY object_id"
    params = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query})
    url = ENDPOINT + "?" + params
    key = hashlib.sha256(query.encode()).hexdigest()[:16]
    snapshot = directory / (key + ".csv")
    for attempt in range(3):
        try:
            if snapshot.exists():
                content = snapshot.read_bytes()
            else:
                request = Request(url, headers={"User-Agent": "Open-Space-Euclid-Q1/1.0"})
                with urlopen(request, timeout=110) as stream:
                    content = stream.read()
            reader = csv.DictReader(StringIO(content.decode()))
            if set(reader.fieldnames or []) != set(FIELDS):
                raise ValueError("PHZ response has unexpected columns")
            rows = list(reader)
            if len({r["object_id"] for r in rows}) != len(rows) or any(r["object_id"] not in ids for r in rows):
                raise ValueError("PHZ returned duplicate or unexpected object_id")
            snapshot.write_bytes(content)
            return rows, {"query": query, "url": url, "requested": len(ids), "rows": len(rows),
                          "sha256": hashlib.sha256(content).hexdigest()}
        except Exception:
            if attempt == 2:
                raise
            sleep(2 ** attempt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    with (args.input_dir / "q1_discovery_engine_lens_catalog.csv").open() as stream:
        ids = {r["object_id"] for r in csv.DictReader(stream)
               if r["subset"] == "discovery_engine" and r["grade"] in {"A", "B"}}
    with (args.input_dir / "agilelens.tsv").open() as stream:
        lines = [line for line in stream if line.strip() and not line.startswith("#")]
    for row in csv.DictReader([lines[0], *lines[3:]], delimiter="\t"):
        ids.add(str(int(row["cutoutName"].strip().split("_")[-1])))
    ids = sorted(ids, key=int)
    cache = args.input_dir / "phz_batches"
    cache.mkdir(exist_ok=True)
    batches = [ids[i:i+55] for i in range(0, len(ids), 55)]
    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(lambda batch: ajax(batch, cache), batches))
    rows = [r for found, _ in results for r in found]
    with (args.output / "phz_merged.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: int(r["object_id"])))
    (args.output / "phz_queries.json").write_text(json.dumps([item for _, item in results], indent=2) + "\n")
    print(f"{len(rows)} returned for {len(ids)} exact IDs in {len(batches)} queries")


if __name__ == "__main__":
    main()
