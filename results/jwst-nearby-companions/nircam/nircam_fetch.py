"""Fetch only the declared public files, verify SHA-256, and publish atomically."""

import argparse
import csv
import hashlib
import tempfile
from pathlib import Path
from urllib.request import urlopen


def digest(path):
    with path.open("rb") as data:
        return hashlib.file_digest(data, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=Path("files"))
    parser.add_argument(
        "--manifest", type=Path, default=Path(__file__).with_name("mcc124_inputs.csv")
    )
    args = parser.parse_args()
    args.inputs.mkdir(parents=True, exist_ok=True)
    manifest = args.manifest
    with manifest.open() as stream:
        for row in csv.DictReader(stream):
            target = args.inputs / row["filename"]
            if target.exists() and digest(target) == row["sha256"]:
                continue
            with tempfile.NamedTemporaryFile(dir=args.inputs, delete=False) as tmp:
                staged = Path(tmp.name)
            try:
                with (
                    urlopen(row["url"], timeout=180) as result,
                    staged.open("wb") as output,
                ):
                    while block := result.read(1024 * 1024):
                        output.write(block)
                if digest(staged) != row["sha256"]:
                    raise ValueError(
                        f"{target.name} differs from the recorded revision; inspect archive reprocessing"
                    )
                staged.replace(target)
            finally:
                staged.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
