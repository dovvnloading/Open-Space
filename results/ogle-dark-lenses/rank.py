"""Rank screened events by the lens mass implied at a reference proper motion.

Usage: python rank.py SCREEN.csv OUT.csv

Selection: converged fit, 50 < tE < 1500 d, sigma(tE)/tE < 0.25, both parallax
component uncertainties < 0.1, |u0| < 1 (peak magnification above 1.34), positive
source flux and a blend fraction above -0.2. The rank
statistic is the lens mass for a relative proper motion of 4 mas/yr when the
parallax is raised by twice its larger component uncertainty:

    M_ref = 4 mas/yr * tE / (kappa * (|piE| + 2 sigma)),  kappa = 8.144 mas/Msun.

It is a conservative ordering statistic, not a mass estimate.
"""
import sys

import numpy as np
import pandas as pd

KAPPA = 8.144


def ranked(path, te_precision=0.25):
    d = pd.read_csv(path)
    counts = {"events": len(d), "fitted": int((d.status == "ok").sum())}
    d = d[d.status == "ok"].copy()
    d["pie"] = np.hypot(d.pi_n, d.pi_e)
    d["sigma_pie"] = d[["e_pi_n", "e_pi_e"]].max(axis=1)
    d["blend_fraction"] = d.fb / (d.fs + d.fb)
    d["m_ref"] = 4 * d.te / 365.25 / KAPPA / (d.pie + 2 * d.sigma_pie)
    keep = ((d.te > 50) & (d.te < 1500) & (d.e_te / d.te < te_precision) & (d.sigma_pie < 0.1)
            & (d.u0.abs() < 1.0) & (d.fs > 0) & (d.blend_fraction > -0.2))
    counts["pass_quality"] = int(keep.sum())
    return d[keep].sort_values("m_ref", ascending=False), counts, d


def main():
    table, counts, _ = ranked(sys.argv[1])
    columns = ["event", "ra", "dec", "n", "t0", "te_0", "te", "e_te", "u0", "pi_n", "e_pi_n", "pi_e", "e_pi_e",
               "pie", "dchi2", "blend_fraction", "i_base", "last_hjd", "m_ref"]
    table[columns].to_csv(sys.argv[2], index=False, float_format="%.5g")
    print(counts)
    print(table[columns].head(25).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
