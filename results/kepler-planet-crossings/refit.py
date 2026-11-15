"""Reoptimize all per-event nuisance parameters for four projected rotations."""

import pandas as pd
from replay import HERE, GEOMETRIES, build_events


def main():
    rows = []
    for column in ("pdc", "sap"):
        for angle in (58, 60, 62, 64):
            events = build_events(GEOMETRIES["ShortCadence"], column, "best", nonlinear=(angle, 1))
            for epoch, delta, null, alt in events:
                rows.append([column, angle, epoch, delta, *null[:2], *alt[:2]])
    columns = ["col", "ang", "ep", "delta", "null_c", "null_d", "alt_c", "alt_d"]
    pd.DataFrame(rows, columns=columns).to_csv(HERE / "nonlinear.csv", index=False)


if __name__ == "__main__":
    main()
