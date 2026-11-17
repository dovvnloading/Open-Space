"""Fetch native B-camera barycentric factors for completed spectral shards."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import json
from pathlib import Path
from frame_headers import get_header


def headers(scan, output, cache):
    with gzip.open(scan, "rt") as f:
        keys = {(e["night"], e["expid"], "b" + str(e["petal_loc"]))
                for line in f for e in json.loads(line)["meta"]["exposures"]}
    cache.mkdir(exist_ok=True,parents=True)
    def cached(key):
        path = cache / f"{key[0]}-{key[1]}-{key[2]}.json"
        if path.exists():
            return json.loads(path.read_text())
        row = get_header(key)
        path.write_text(json.dumps(row))
        return row
    with ThreadPoolExecutor(max_workers=8) as pool:
        rows = list(pool.map(cached, sorted(keys)))
    if not rows:
        return
    partial = output.with_suffix(".partial")
    with partial.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(output)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder",type=Path)
    a=p.parse_args()
    for c in sorted(a.folder.glob("scale-count-*.json")):
        tag=c.stem.removeprefix("scale-count-")
        dst=a.folder / f"scale-headers-{tag}.csv"
        if not dst.exists():
            headers(a.folder / f"scale-scan-{tag}.jsonl.gz",dst,a.folder/"header-prefix-cache")
