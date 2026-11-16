"""Attach an observed selection outcome to eight coordinate control regions."""

import argparse
import csv
from pathlib import Path
import numpy as np
from astropy.io import fits
from simbad_screen import separation


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("catalogue", type=Path)
    ap.add_argument("controls", type=Path)
    ap.add_argument("output_directory", type=Path)
    a = ap.parse_args()
    selected = {}
    for name in ("measured", "boundary", "rescue", "bright"):
        with (a.output_directory / f"{name}.csv").open() as f:
            selected.update({r["SRCID"]: name for r in csv.DictReader(f)})
    with a.controls.open() as f:
        controls = list(csv.DictReader(f))
    out = []
    with fits.open(a.catalogue, memmap=True) as hdus:
        d = hdus[1].data
        ras = d["RA"].astype(float)
        decs = d["DEC"].astype(float)
        for row in controls:
            rr, dd = float(row["ra"]), float(row["dec"])
            pre = np.abs(decs - dd) < 60.0 / 3600.0
            nearby = [
                (separation(rr, dd, ras[i], decs[i]), i) for i in np.flatnonzero(pre)
            ]
            radius, i = min(nearby)
            if radius > 60:
                raise ValueError(f"No 5XMM source within control search radius: {row}")
            src = d[i]
            out.append(
                dict(
                    control=row["id"],
                    SIMBAD_RA=rr,
                    SIMBAD_DEC=dd,
                    SRCID=int(src["SRCID"]),
                    IAUNAME=src["IAUNAME"],
                    distance_arcsec=round(radius, 4),
                    branch=selected.get(str(src["SRCID"]), "outside"),
                    SUM_FLAG=src["SUM_FLAG"],
                    EXTENT=src["EXTENT"],
                    DIST_NN=src["DIST_NN"],
                    STACK_DET_ML=src["STACK_DET_ML"],
                    EP_CTS=src["EP_CTS"],
                    RADEC_ERR=src["RADEC_ERR"],
                    HR1=src["EP_HR1"],
                    HR1_ERR=src["EP_HR1_ERR"],
                    HR2=src["EP_HR2"],
                    HR2_ERR=src["EP_HR2_ERR"],
                    VAR_FLAG=bool(src["VAR_FLAG"]),
                    VAR_PROB=src["VAR_PROB"],
                    FRATIO=src["FRATIO"],
                    FRATIO_ERR=src["FRATIO_ERR"],
                    B4_SN=src["EP_4_FLUX"] / src["EP_4_FLUX_ERR"],
                    B5_SN=src["EP_5_FLUX"] / src["EP_5_FLUX_ERR"],
                )
            )
    with (a.output_directory / "control_recovery.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, list(out[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    main()
