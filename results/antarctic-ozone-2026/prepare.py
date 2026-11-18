"""Build a lossless polar subset of the pinned NASA ozone maps."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

import numpy as np

from source import read_grid

ROOT = Path(__file__).resolve().parent


def run(data, output, metadata):
    source = json.loads((ROOT / "sources.json").read_bytes())
    if source["unavailable"]:
        raise ValueError("The selected date and sensor manifest is incomplete")
    saved = {}
    with ZipFile(output, "w") as bundle:
        for entry in source["records"]:
            path = data / entry["file"]
            body = path.read_bytes()
            if len(body) != entry["bytes"] or sha256(body).hexdigest() != entry["raw_sha256"]:
                raise ValueError(f"Source does not match the receipt: {path}")
            _, grid = read_grid(path, expected_date=entry["date"], expected_sensor=entry["sensor"].upper())
            member = f"{entry['date']}-{entry['sensor']}.u16"
            cap = np.asarray(grid[:50], dtype="<u2")
            bundle.writestr(ZipInfo(member), cap.tobytes(), compress_type=ZIP_DEFLATED, compresslevel=9)
            saved[member] = {
                "source_file": path.name,
                "first_header_line": body.splitlines()[0].decode("ascii"),
                "polar_sha256": sha256(cap.tobytes()).hexdigest(),
            }
    metadata.write_text(json.dumps({
        "units": "Dobson units",
        "shape": [50, 360],
        "representation": "row-major little-endian unsigned 16-bit; south-to-north, west-to-east",
        "values": "exact native positive values and zeros retained",
        "members": saved,
    }, indent=2) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    args = parser.parse_args()
    run(args.data_root, args.output, args.metadata)
