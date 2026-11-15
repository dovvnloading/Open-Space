"""Read the NASA K2C9 tarball and verify the four original tables and selected rows."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import tarfile
from urllib.request import urlopen

URL = "https://exoplanet.ipac.caltech.edu/data/Contributed/K2C9/K2C9.tar.gz"
SHA = "69ebcfdd40d85a1c911e314b359fe2f21444da0deac64bbe98d7d48c21a0d2eb"
TABLES = {
    "exofop_k2c9_events_20220315_135352.csv": "bddbdafe0b5f0d6395de3a1abe3c17ae7febd86a8a806dd372649353965c3c1b",
    "k2c9_event_param_2022-03-15_13:25:32_nulls_corrected.txt": "ab2e297036a0f32a358bd2769e88c16b48ab7ab2019cc0b151948e8ad58b2da9",
    "k2c9_files_2022-03-15_13:25:32.txt": "80ca1fdcaf4ce29010b7a926ba9be0926b32dc7e90d6e46670bc928b0e16b9d5",
    "k2c9_changelog_2022-03-15_13:25:32.txt": "799683ee3f357d9c5705ce7183805354c1c9a3438e641dcc861578b56826f5c7",
}
EID = "K2C9-R-0367"
OGLE = "OGLE-2016-BLG-0935"


class HashingReader:
    def __init__(self, stream):
        self.stream = stream
        self.digest = hashlib.sha256()

    def read(self, n):
        chunk = self.stream.read(n)
        self.digest.update(chunk)
        return chunk


def verified_rows(tables):
    out = {}
    if set(tables) != set(TABLES):
        raise ValueError("Missing original table")
    for name in sorted(tables):
        contents = tables[name]
        if hashlib.sha256(contents).hexdigest() != TABLES[name]:
            raise ValueError(f"Unexpected table hash: {name}")
        lines = io.StringIO(contents.decode("utf-8-sig"))
        reader = csv.DictReader(lines, delimiter="," if name.endswith(".csv") else "|")
        matches = [row for row in reader if any(
            EID in (v or "") or OGLE in (v or "") for v in row.values()
        )]
        out[name] = matches
    return out


def extract(stream):
    source = HashingReader(stream)
    tables = {}
    with tarfile.open(fileobj=source, mode="r|gz") as archive:
        for member in archive:
            name = member.name.removeprefix("./")
            if member.isfile() and name in TABLES:
                tables[name] = archive.extractfile(member).read()
    while source.read(1 << 20):
        pass
    if source.digest.hexdigest() != SHA:
        raise ValueError("Unexpected archive SHA-256")
    return verified_rows(tables)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="optional local K2C9.tar.gz")
    args = parser.parse_args()
    with (args.archive.open("rb") if args.archive else urlopen(URL, timeout=150)) as stream:
        rows = extract(stream)
    print(json.dumps(rows, indent=2, ensure_ascii=False))
