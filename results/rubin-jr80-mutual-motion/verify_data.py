"""Verify the input snapshot; optionally repeat public TAP queries to a new directory."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from cutouts import DATA

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh-dir", type=Path)
    args = parser.parse_args()
    meta = json.loads((DATA/"provenance.json").read_text())
    for file, digest in meta["frozen_sha256"].items():
        if hashlib.sha256((DATA/file).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Unexpected bytes in {file}")
    if args.refresh_dir:
        if args.refresh_dir.resolve() == DATA.resolve():
            raise ValueError("Use another directory for live responses")
        args.refresh_dir.mkdir(parents=True, exist_ok=True)
        for name in ("screen","orbits","target","peer_pool"):
            query = (DATA/f"{name}.adql").read_text().strip()
            url = meta["tap_sync"]+"?"+urlencode({"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":query})
            with urlopen(url, timeout=180) as response:
                result=response.read()
            header = (DATA/f"{name}.csv").read_bytes().splitlines()[0]
            if result.splitlines()[0] != header:
                raise ValueError(f"Unexpected TAP response for {name}")
            (args.refresh_dir/f"{name}.csv").write_bytes(result)
    print("Verified frozen tables, queries and image manifests.")

if __name__ == "__main__":
    main()
