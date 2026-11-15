"""Collapse NEOWISE-R single exposures into one flux per sky visit and source.

Usage: python epochs.py BATCHDIR OUT.csv [--radius ARCSEC]

A visit is a group of exposures separated from the next by more than 60 days.
Exposures need qual_frame > 0, qi_fact > 0, saa_sep > 5 deg and an unmasked
Moon flag for the band. The visit flux is the median of the exposure fluxes;
its uncertainty is 1.4826 MAD / sqrt(n), with the median catalog uncertainty
/ sqrt(n) as a floor. Fluxes are in mJy (Vega zero points 309.54 and 171.787 Jy).
"""
import argparse
import csv
import gzip
import multiprocessing as mp
import pathlib

import numpy as np
import pandas as pd

NAMES = ["cntr", "dist", "pang", "ra_in", "dec_in", "ra", "dec", "clon", "clat", "mjd", "w1", "e1", "w2", "e2",
         "chi1", "chi2", "qual_frame", "qi_fact", "saa_sep", "moon_masked", "cc_flags", "sso_flg"]
ZERO = {"1": 309.54e3, "2": 171.787e3}


def read_batch(path):
    with gzip.open(path, "rt") as handle:
        skip = 0
        for line in handle:
            if line.startswith("\\") or line.startswith("|"):
                skip += 1
            else:
                break
    return pd.read_csv(path, sep=r"\s+", skiprows=skip, names=NAMES, na_values=["null"],
                       dtype={"moon_masked": str, "cc_flags": str}, compression="gzip")


def visits(group, radius):
    group = group[(group.dist <= radius) & (group.qual_frame > 0) & (group.qi_fact > 0) & (group.saa_sep > 5)]
    if len(group) == 0:
        return []
    group = group.sort_values("mjd")
    label = np.cumsum(np.diff(group.mjd.values, prepend=group.mjd.values[0]) > 60)
    rows = []
    for _, visit in group.groupby(label):
        row = {"mjd": float(visit.mjd.median())}
        moon = visit.moon_masked.fillna("00").str.zfill(2)
        for band, index in (("1", 0), ("2", 1)):
            use = visit[f"w{band}"].notna() & visit[f"e{band}"].notna() & (moon.str[index] == "0")
            mags = visit.loc[use, f"w{band}"].values
            n = len(mags)
            if n < 4:
                row.update({f"f{band}": np.nan, f"s{band}": np.nan, f"n{band}": n})
                continue
            flux = ZERO[band] * 10 ** (-0.4 * mags)
            median = np.median(flux)
            scatter = 1.4826 * np.median(np.abs(flux - median))
            catalog = median * 0.921 * np.median(visit.loc[use, f"e{band}"].values)
            row.update({f"f{band}": median, f"s{band}": max(scatter, catalog) / np.sqrt(n), f"n{band}": n})
        row["dist"] = float(visit.dist.median())
        rows.append(row)
    return rows


def one_batch(job):
    path, radius = job
    ids = (path.parent / path.name.replace(".tbl.gz", ".ids")).read_text().splitlines()
    out = []
    for cntr, group in read_batch(path).groupby("cntr"):
        for row in visits(group, radius):
            out.append([ids[int(cntr) - 1], f"{row['mjd']:.2f}"] +
                       [f"{row[k]:.5g}" for k in ("f1", "s1", "n1", "f2", "s2", "n2", "dist")])
    return path.name, out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("batchdir")
    parser.add_argument("out")
    parser.add_argument("--radius", type=float, default=3.0)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    jobs = [(path, args.radius) for path in sorted(pathlib.Path(args.batchdir).glob("batch_*.tbl.gz"))]
    with mp.Pool(args.workers) as pool, open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", "mjd", "f1", "s1", "n1", "f2", "s2", "n2", "dist"])
        for name, rows in pool.imap(one_batch, jobs):
            writer.writerows(rows)
            print(name, len(rows), flush=True)


if __name__ == "__main__":
    main()
