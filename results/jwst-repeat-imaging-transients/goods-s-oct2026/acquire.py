"""Acquire the frozen GOODS-S catalogue snapshot without changing its receipts."""
import argparse
import csv
import hashlib
from pathlib import Path, PurePosixPath
import re
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
BASE = "https://stpubdata.s3.amazonaws.com/"


def receipts(pid, directory=HERE):
    if str(pid) not in {"8060", "12577"}:
        raise ValueError("Unsupported frozen program")
    prefix = f"jwst/public/jw{int(pid):05d}/L3/t/"
    with (directory / f"{pid}_index.csv").open(newline="") as f:
        index = list(csv.DictReader(f))
    with (directory / f"{pid}_manifest.csv").open(newline="") as f:
        manifest = list(csv.DictReader(f))
    indexed = {}
    for row in index:
        key = row["key"]
        if key in indexed or not key.startswith(prefix) or ".." in PurePosixPath(key).parts or "\\" in key:
            raise ValueError("Invalid or duplicate index key")
        indexed[key] = int(row["bytes"])
    cats = {k: v for k, v in indexed.items() if k.endswith("_cat.ecsv")}
    found, paths = {}, set()
    for row in manifest:
        url, local, sha = row["url"], row["local_path"], row["sha256"]
        if not url.startswith(BASE):
            raise ValueError("Invalid source URL")
        key = url[len(BASE):]
        path = PurePosixPath(local)
        if (len(path.parts) != 2 or path.parts[0] != str(pid) or path.name != PurePosixPath(key).name
                or "\\" in local or ".." in path.parts or local != path.as_posix()):
            raise ValueError("Invalid local path")
        size = int(row["bytes"])
        if key in found or local in paths or key not in cats or cats[key] != size or size <= 0:
            raise ValueError("Duplicate or inconsistent manifest member")
        if not re.fullmatch(r"[0-9a-f]{64}", sha):
            raise ValueError("Invalid SHA-256")
        found[key] = row
        paths.add(local)
    if not cats or set(found) != set(cats):
        raise ValueError("Manifest does not cover frozen catalog index")
    return list(found.values())


def acquire(pid, cache_root, directory=HERE, fetch=None):
    rows = receipts(pid, directory)
    cache_root = Path(cache_root)
    if fetch is None:
        def fetch(url):
            with urlopen(url, timeout=90) as stream:
                return stream.read()
    pending = []
    for row in rows:
        path = cache_root / Path(*PurePosixPath(row["local_path"]).parts)
        data = path.read_bytes() if path.exists() else fetch(row["url"])
        if len(data) != int(row["bytes"]) or hashlib.sha256(data).hexdigest() != row["sha256"]:
            raise ValueError(f"Frozen catalogue mismatch: {row['local_path']}")
        if not path.exists():
            pending.append((path, data))
    # No cache product or receipt is written until every input has passed.
    for path, data in pending:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    return len(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("pid", choices=("8060", "12577"))
    p.add_argument("--cache-root", type=Path, default=Path.cwd(), help="Root containing PID directories (use D: for downloads)")
    args = p.parse_args()
    print(f"Verified {acquire(args.pid, args.cache_root)} frozen catalogues; receipts unchanged")


if __name__ == "__main__":
    main()
