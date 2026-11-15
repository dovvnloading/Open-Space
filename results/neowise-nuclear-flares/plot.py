"""Plot NEOWISE visit light curves for a list of sources.

Usage: python plot.py EPOCHS.csv OUT.png ID [ID ...] [--titles "A|B|..."]
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("epochs")
    parser.add_argument("out")
    parser.add_argument("ids", nargs="+")
    parser.add_argument("--titles", default="")
    args = parser.parse_args()
    table = pd.read_csv(args.epochs, dtype={"id": str})
    titles = args.titles.split("|") if args.titles else args.ids
    columns = min(3, len(args.ids))
    rows = int(np.ceil(len(args.ids) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(5.2 * columns, 3.4 * rows), squeeze=False)
    for axis, name, title in zip(axes.flat, args.ids, titles):
        group = table[table.id == name].sort_values("mjd")
        year = 2000 + (group.mjd - 51544.5) / 365.25
        axis.errorbar(year, group.f1, yerr=group.s1, fmt="o-", ms=4, color="C0", label="W1 (3.4 um)")
        axis.errorbar(year, group.f2, yerr=group.s2, fmt="s-", ms=4, color="C3", label="W2 (4.6 um)")
        axis.set_title(title, fontsize=10)
        axis.set_ylabel("flux (mJy)")
        axis.set_xlabel("year")
        axis.set_ylim(bottom=0)
        axis.legend(fontsize=7, loc="upper left")
    for axis in axes.flat[len(args.ids):]:
        axis.set_visible(False)
    fig.tight_layout()
    fig.savefig(args.out, dpi=110)


if __name__ == "__main__":
    main()
