"""Verify the actual resampled image arrays independently of NPZ container metadata."""
import csv
from hashlib import sha256
from pathlib import Path
import sys

import numpy as np


if __name__ == "__main__":
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("stamps")
    with Path(__file__).with_name("pixel-digests.csv").open() as manifest:
        for row in csv.DictReader(manifest):
            with np.load(folder / row["file"]) as stamp:
                data = np.asarray(stamp["image"], dtype="<f4", order="C")
            if sha256(data.tobytes()).hexdigest() != row["sha256"]:
                raise ValueError(f"Pixel digest mismatch in {row['file']}")
    print("All resampled pixel digests verified")
