"""Collect per-integration CSVs after running the numerical experiments."""

from pathlib import Path
import csv

WORK = Path(__file__).resolve().parent / "work"


def collect(pattern, name, expected):
    paths = sorted(WORK.glob(pattern))
    if len(paths) != expected:
        raise RuntimeError(f"{pattern} matched {len(paths)} runs, expected {expected}")
    with (WORK / name).open("w") as stream:
        writer = csv.writer(stream)
        for index, path in enumerate(paths):
            with path.open() as source:
                reader = csv.reader(source)
                header = next(reader)
                if index == 0:
                    writer.writerow(header)
                writer.writerows(reader)


if __name__ == "__main__":
    collect("current_*.csv", "mpc-clones.csv", 20)
    collect("deep_*.csv", "vz98-deep.csv", 6)
    collect("target_*.csv", "old-targets.csv", 12)
