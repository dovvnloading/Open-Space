"""Nominate sources with no same-filter previous counterpart from independent long bands."""

import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from catalog_screen import readcat, project, match

r = []
for obs in [1, 2]:
    for curr, prev in [(12577, 8060), (8060, 12577)]:
        bits = []
        for band in ["277", "356"]:
            c = readcat(curr, obs, band)
            a = readcat(prev, obs, band)
            xy = project(c)
            tree = cKDTree(project(a))
            d, j = tree.query(xy, k=4)
            cut = (
                (d[:, 0] > 0.3)
                & (d[:, -1] < 7)
                & (c.aper_total_flux / c.aper_total_flux_err >= 12)
                & ~c.is_extended
            )
            bits.append(c[cut].copy().reset_index(drop=True))
        if any(len(t) == 0 for t in bits):
            continue
        u, j = match(project(bits[0]), project(bits[1]), 0.10)
        for i in np.flatnonzero(u):
            t = bits[0].iloc[i]
            s = bits[1].iloc[j[i]]
            r.append(
                dict(
                    obs=obs,
                    present=curr,
                    ra=(t["sky_centroid.ra"] + s["sky_centroid.ra"]) / 2,
                    dec=(t["sky_centroid.dec"] + s["sky_centroid.dec"]) / 2,
                    f277_njy=t.aper_total_flux * 1e9,
                    f356_njy=s.aper_total_flux * 1e9,
                    f277_snr=t.aper_total_flux / t.aper_total_flux_err,
                    f356_snr=s.aper_total_flux / s.aper_total_flux_err,
                )
            )
out = pd.DataFrame(r)
out["site"] = 1000 + np.arange(len(out))
out.to_csv("onsets.csv", index=False)
print(out.round(3).to_string(index=False))
