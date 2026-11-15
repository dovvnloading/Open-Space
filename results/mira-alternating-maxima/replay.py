"""Rank every star by its longest strictly alternating run of consecutive maxima."""

import csv
import gzip
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load_rows(filename):
    with (HERE / filename).open(newline="") as stream:
        return list(csv.DictReader(stream))


def load_series():
    series = defaultdict(dict)
    with gzip.open(HERE / "cycle-magnitudes.csv.gz", "rt") as stream:
        for row in csv.DictReader(stream):
            name, cycle = row["star"], int(row["cycle"])
            if cycle in series[name]:
                raise ValueError(f"Duplicate cycle {name} {cycle}")
            series[name][cycle] = int(row["mag_tenths"])
    return series


def segments(series):
    segment = []
    for cycle, mag in sorted(series.items()):
        if segment and cycle != segment[-1][0] + 1:
            yield segment
            segment = []
        segment.append((cycle, mag))
    if segment:
        yield segment


def longest(segment):
    start = best_start = best_end = 0
    previous = 0
    for index in range(1, len(segment)):
        change = segment[index][1] - segment[index - 1][1]
        sign = (change > 0) - (change < 0)
        if sign == 0:
            start = index
        elif sign == previous:
            start = index - 1
        if index - start > best_end - best_start:
            best_start, best_end = start, index
        previous = sign
    return best_end - best_start + 1, segment[best_start][0], segment[best_end][0]


if __name__ == "__main__":
    import sys

    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(("star", "n_peaks", "longest_run", "first_cycle", "last_cycle"))
    for star, series in sorted(load_series().items()):
        winner = max(longest(segment) for segment in segments(series))
        writer.writerow([star, len(series), *winner])
