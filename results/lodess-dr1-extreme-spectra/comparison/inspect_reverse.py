"""Attach the survey photometry to each VLSSr-parent catalogue nonmatch."""

import argparse
from pathlib import Path

import pandas as pd


def join_results(root):
    out = root / "results"
    phot = pd.read_csv(out / "vlssr_nondetections.csv")
    case = (
        pd.read_csv(out / "reverse_positions.csv")
        .merge(phot, on=["name", "ra", "dec"])
        .set_index("row")
    )
    for name in ["gmrtas150m", "nvss", "rrs8c38mhz", "wenss", "gleamegcat", "lotss"]:
        path = root / "reverse_results" / f"match_{name}.csv.gz"
        case = case.join(pd.read_csv(path, index_col="row"), validate="many_to_one")
    case.to_csv(out / "reverse_counterparts.csv.gz", index=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("comparison", type=Path)
    args = p.parse_args()
    join_results(args.comparison)
