"""Read a spectral manifest or measurements once per target and native image."""
import csv
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit


def read_frames(paths):
    rows = []
    seen = set()
    for path in paths:
        with Path(path).open() as stream:
            for row in csv.DictReader(stream):
                product = PurePosixPath(urlsplit(row["cutout_url"]).path).name
                key = row["target"], product
                if key in seen:
                    raise ValueError(f"Repeated native image for target {key[0]}: {product}")
                seen.add(key)
                rows.append(row)
    return rows
