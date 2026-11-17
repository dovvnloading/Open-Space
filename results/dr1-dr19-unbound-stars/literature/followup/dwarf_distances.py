"""Evaluate single-dwarf photometric distances, conditional on the adopted luminosity class."""
import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def ps1_to_sdss(g, r, i):
    # Tonry et al. 2012, Table 6, quadratic transformations as functions of observed PS1 g-r.
    x = g - r
    return (g + 0.013 + 0.145*x + 0.019*x*x,
            r - 0.001 + 0.004*x + 0.007*x*x,
            i - 0.005 + 0.011*x + 0.010*x*x)


def main():
    rows = []
    with (HERE/"photometric_inputs.csv").open() as f:
        for r in csv.DictReader(f):
            g, red, i = map(float, (r["g"],r["r"],r["i"]))
            if r["survey"] == "PS1":
                g, red, i = ps1_to_sdss(g,red,i)
            g -= float(r["ag"])
            red -= float(r["ar"])
            i -= float(r["ai"])
            x = g - i
            # Ivezić et al. 2008, Eqs A2 and A7. The metal-poor calibrator M15 is near -2.5 dex.
            m = max(float(r["feh_desi"]), -2.5)
            mr = -0.56 + 14.32*x - 12.97*x**2 + 6.127*x**3 - 1.267*x**4 + 0.0967*x**5 - 1.11*m - 0.18*m*m
            if not (0.2 < x < 4.0):
                continue
            d = 10**((red-mr-10)/5)
            rows.append([r["source_id"],r["survey"],r["object_id"],r["pos_offset_arcsec"],f"{x:.5f}",
                         f"{red:.5f}",f"{m:.5f}",f"{mr:.5f}",f"{d:.5f}",f"{d*math.sqrt(2):.5f}"])
    with (HERE/"dwarf_distances.csv").open("w",newline="") as f:
        w=csv.writer(f,lineterminator="\n")
        w.writerow(["source_id","survey","object_id","pos_offset_arcsec","g_i_dered","r_dered",
                    "feh_applied","Mr_single","distance_single_kpc","distance_equal_luminosity_pair_kpc"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
