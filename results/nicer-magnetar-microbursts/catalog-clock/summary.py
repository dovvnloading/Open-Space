"""Summarize paired windows in the published NICER burst list, offline."""

import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent


def summarize(rows, width):
    literal = [int(r[f"{width}_literal"]) for r in rows]
    shifted = [int(r[f"{width}_plus69"]) for r in rows]
    wins = sum(a > b for a, b in zip(shifted, literal))
    ties = sum(a == b for a, b in zip(shifted, literal))
    return (len(rows), wins, ties, len(rows)-wins-ties, sum(literal), sum(shifted),
            median([a-b for a, b in zip(shifted, literal)]))


def main():
    with (ROOT / "windows.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    print(len(rows), "catalog rows in", len({r["obsid"] for r in rows}), "observations")
    for width in ("narrow", "wide"):
        selected = [r for r in rows if r[f"{width}_literal"] and r[f"{width}_plus69"]]
        print(width, "pairs, wins, ties, losses, literal hits, shifted hits, median difference")
        print("all", summarize(selected, width))
        group = defaultdict(list)
        obs = defaultdict(list)
        for r in selected:
            group[r["source"]].append(r)
            obs[r["obsid"]].append(r)
        for name, sample in sorted(group.items()):
            print(name, summarize(sample, width))
        strong = sum(summarize(v, width)[-1] > 0 for v in obs.values())
        print("observations with positive median difference", strong, "of", len(obs))
        four = [r for r in selected if r[f"{width}_minus69"] and r[f"{width}_plus138"]]
        print("four-way complete",len(four),{d:sum(int(r[f"{width}_{d}"]) for r in four)
                                               for d in ("minus69","literal","plus69","plus138")})
        print("single-source storm excluded", summarize([r for r in selected if r["source"] != "SGR 1935+2154"],width))


if __name__ == "__main__":
    main()
