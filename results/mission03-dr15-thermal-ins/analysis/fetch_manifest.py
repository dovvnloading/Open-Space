"""Fetch and verify the archive inputs used by the pointed-observation audit."""

import argparse
import csv
import hashlib
from pathlib import Path
from urllib.request import urlopen
from tempfile import NamedTemporaryFile
from concurrent.futures import ThreadPoolExecutor


def restore(row, directory):
    path = directory / row["file"]

    def matches(p):
        with p.open("rb") as f:
            return hashlib.file_digest(f, "sha256").hexdigest() == row["sha256"]

    if path.exists() and matches(path):
        return
    with NamedTemporaryFile(dir=directory, delete=False) as tmp:
        tmp_path = Path(tmp.name)
        try:
            with urlopen(row["url"], timeout=120) as r:
                while blk := r.read(1024 * 1024):
                    tmp.write(blk)
            tmp.flush()
            if not matches(tmp_path):
                raise ValueError(f"Hash mismatch for {path.name}")
            tmp_path.replace(path)
        finally:
            tmp_path.unlink(missing_ok=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument(
        "--group",
        default="J085618",
        help="Source tag, control name or erass1/hi4pi/mekal",
    )
    a = p.parse_args()
    a.directory.mkdir(parents=True, exist_ok=True)
    with Path(__file__).with_name("archive_manifest.csv").open() as f:
        rows = [r for r in csv.DictReader(f) if r["group"] == a.group]
    if not rows:
        raise ValueError("No such manifest group")
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(lambda r: restore(r, a.directory), rows))


if __name__ == "__main__":
    main()
