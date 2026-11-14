"""Extract the measurement columns from the original mission products."""
import argparse
import csv
import hashlib
from pathlib import Path

from astropy.io import fits

ROOT = Path(__file__).resolve().parent
FIELDS = ("TSTART", "TSTOP", "BARYTIME", "FREQUENCY", "FREQUENCY_ERR",
          "FMIN", "FMAX", "NPROFILES", "DETECTED", "AMPLITUDE", "AMPSIG")


def extract(gbm, maxi, output):
    with fits.open(gbm) as hdus:
        history = hdus["PSR_HIST"].data
        with (output / "gbm.csv").open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow([key.lower() for key in FIELDS])
            for row in history:
                writer.writerow([row[key].item() for key in FIELDS])
    with maxi.open() as source, (output / "maxi.csv").open("w", newline="") as dest:
        writer = csv.writer(dest, lineterminator="\n")
        writer.writerow(["mjd", "flux_2_20", "error_2_20", "flux_2_4", "error_2_4",
                         "flux_4_10", "error_4_10", "flux_10_20", "error_10_20"])
        for line in source:
            fields = line.split()
            if fields and 61046 <= float(fields[0]) < 61241:
                writer.writerow(fields)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gbm", type=Path)
    parser.add_argument("maxi", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for path in [args.gbm, args.maxi]:
        print(hashlib.sha256(path.read_bytes()).hexdigest(), path)
    extract(args.gbm, args.maxi, args.output)
