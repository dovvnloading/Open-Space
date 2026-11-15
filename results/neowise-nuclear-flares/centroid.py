"""Measure whether the infrared centroid moved during a flare.

Usage: python centroid.py BATCHDIR FLARES.csv OUT.csv

For each flare the single-exposure profile-fit positions are split at the
onset date. A transient at the nucleus leaves the centroid in place; one in
the disc pulls the flare-epoch positions towards itself. The shift of the
blended centroid is multiplied by (total flux / excess flux) at the peak to
estimate the offset of the transient itself.
"""
import csv
import pathlib
import sys

import numpy as np
import pandas as pd

import epochs


def main():
    batchdir = pathlib.Path(sys.argv[1])
    flares = pd.read_csv(sys.argv[2], dtype={"id": str}).set_index("id")
    where = {}
    for path in sorted(batchdir.glob("batch_*.ids")):
        for i, name in enumerate(path.read_text().splitlines(), 1):
            if name in flares.index:
                where.setdefault(path, {})[i] = name
    rows = []
    for path, wanted in where.items():
        table = epochs.read_batch(path.with_name(path.name.replace(".ids", ".tbl.gz")))
        for cntr, name in wanted.items():
            flare = flares.loc[name]
            g = table[(table.cntr == cntr) & (table.qual_frame > 0) & (table.qi_fact > 0) & (table.saa_sep > 5)
                      & (table.dist <= 3.0) & table.w1.notna()]
            before = g[g.mjd < flare.onset_mjd - 30]
            during = g[(g.mjd > flare.onset_mjd - 30) & (g.mjd < flare.peak_mjd + 400)]
            if len(before) < 20 or len(during) < 8:
                continue
            cosd = np.cos(np.radians(g.dec.median()))
            d_east = (during.ra.median() - before.ra.median()) * cosd * 3600
            d_north = (during.dec.median() - before.dec.median()) * 3600
            err_e = np.hypot(1.4826 * (during.ra - during.ra.median()).abs().median() * cosd * 3600 / np.sqrt(len(during)),
                             1.4826 * (before.ra - before.ra.median()).abs().median() * cosd * 3600 / np.sqrt(len(before)))
            err_n = np.hypot(1.4826 * (during.dec - during.dec.median()).abs().median() * 3600 / np.sqrt(len(during)),
                             1.4826 * (before.dec - before.dec.median()).abs().median() * 3600 / np.sqrt(len(before)))
            excess_share = flare.peak_excess_w1_mjy / (flare.level_w1_mjy + flare.peak_excess_w1_mjy)
            shift = np.hypot(d_east, d_north)
            rows.append({"id": name, "n_before": len(before), "n_during": len(during),
                         "shift_east_arcsec": d_east, "shift_north_arcsec": d_north, "shift_arcsec": shift,
                         "shift_error_arcsec": np.hypot(err_e, err_n),
                         "peak_excess_share_w1": excess_share,
                         "implied_transient_offset_arcsec": shift / max(excess_share, 0.05)})
    with open(sys.argv[3], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in row.items()})
    print(len(rows), "flares measured")


if __name__ == "__main__":
    main()
