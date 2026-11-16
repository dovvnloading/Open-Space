"""Fetch the NASA HEASARC source catalogues and verify the pinned digests."""

import argparse
import csv
import hashlib
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def verify(path, expected):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest() == expected


def acquire(root, row):
    target = root / (row["table"] + ".tdat.gz")
    if target.exists() and verify(target, row["sha256"]):
        return
    with tempfile.NamedTemporaryFile(dir=root, delete=False) as tmp:
        part = Path(tmp.name)
    try:
        subprocess.run(
            [
                "curl",
                "-fsSL",
                "--retry",
                "3",
                "--max-time",
                "900",
                row["url"],
                "-o",
                str(part),
            ],
            check=True,
        )
        if not verify(part, row["sha256"]):
            raise ValueError(f"Digest mismatch for {row['table']}")
        os.replace(part, target)
    finally:
        part.unlink(missing_ok=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("raw_dir", type=Path)
    a = p.parse_args()
    a.raw_dir.mkdir(parents=True, exist_ok=True)
    with (Path(__file__).parent / "archives.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(lambda row: acquire(a.raw_dir, row), rows))
