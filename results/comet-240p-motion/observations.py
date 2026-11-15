"""Read MPC optical positions and match the two components at each observing site."""

from datetime import datetime, timezone

import numpy as np
import pandas as pd


def read_mpc(path):
    rows = []
    for line in path.read_text().splitlines():
        day1 = datetime(int(line[15:19]), int(line[20:22]), 1, tzinfo=timezone.utc)
        jd = 2440587.5 + day1.timestamp() / 86400 + float(line[23:32]) - 1
        h, m, s = map(float, line[32:44].split())
        ra = 15 * (h + m / 60 + s / 3600)
        deg, m, s = map(float, line[45:56].split())
        dec = (deg + m / 60 + s / 3600) * (-1 if line[44] == "-" else 1)
        rows.append((line[11] == "b", line[77:80], int(jd - .5), jd, ra, dec))
    return pd.DataFrame(rows, columns=["secondary", "station", "night", "jd", "ra", "dec"])


def match_components(data):
    """Prefer simultaneous pairs. Otherwise interpolate strictly within the same night."""
    rows = []
    for (site, night), group in data.groupby(["station", "night"]):
        a = group[~group.secondary].sort_values("jd")
        b = group[group.secondary].sort_values("jd")
        common = np.intersect1d(a.jd, b.jd)
        if len(common):
            second = b[b.jd.isin(common)].drop_duplicates("jd")
            kind = "same UTC"
        elif len(a) >= 2 and np.ptp(a.jd) >= .001:
            second = b[b.jd.between(a.jd.min(), a.jd.max())]
            kind = "bracketed"
        else:
            continue
        for obs in second.itertuples():
            ra = np.interp(obs.jd, a.jd, a.ra)
            dec = np.interp(obs.jd, a.jd, a.dec)
            east = 3600 * (obs.ra - ra) * np.cos(np.deg2rad(dec))
            north = 3600 * (obs.dec - dec)
            rows.append((site, night, obs.jd, east, north, kind))
    return pd.DataFrame(rows, columns=["station", "night", "jd", "E", "N", "pairing"])
