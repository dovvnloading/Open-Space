"""Common MPC catalogue and state-vector routines for the 3I–383P calculation."""
from datetime import datetime, timedelta
from pathlib import Path
import re

import numpy as np
from scipy.interpolate import CubicHermiteSpline

INPUTS = Path(__file__).resolve().parent / "inputs"
J2000 = datetime(2000, 1, 1, 12)
AU_KM = 149597870.7
K = 0.01720209895  # Gaussian gravitational constant in au and days


def day(date):
    return (datetime.fromisoformat(date) - J2000).total_seconds() / 86400


def date_string(t):
    return (J2000 + timedelta(days=float(t))).isoformat(timespec="seconds")


def catalog():
    """Read fixed columns of MPC CometEls. Retain distinct input rows."""
    rows = []
    for row in (INPUTS / "CometEls.txt").read_text().splitlines():
        name = row[102:159].strip()
        q, e, arg, node, incl = map(float, row[29:79].split())
        if q > 4.5 or re.search(r"-[A-Z](?:\s|/|$)", name) or "-" in name.split("/")[0]:
            continue
        year, month = int(row[14:18]), int(row[19:21])
        peri = datetime(year, month, 1) + timedelta(days=float(row[22:29]) - 1)
        epoch = datetime.strptime(row[81:89], "%Y%m%d")
        rows.append((name, (peri-J2000).total_seconds()/86400, q, e, arg, node, incl,
                     (epoch-J2000).total_seconds()/86400))
    names = np.array([r[0] for r in rows])
    tp, q, e, arg, node, incl, epochs = np.array([r[1:] for r in rows]).T
    w, o, i = np.deg2rad([arg, node, incl])
    c, s = np.cos, np.sin
    p = np.stack([c(o)*c(w)-s(o)*s(w)*c(i), s(o)*c(w)+c(o)*s(w)*c(i), s(w)*s(i)], axis=1)
    v = np.stack([-c(o)*s(w)-s(o)*c(w)*c(i), -s(o)*s(w)+c(o)*c(w)*c(i), c(w)*s(i)], axis=1)
    return names, tp, q, e, epochs, p, v


def state_spline(a):
    return CubicHermiteSpline(a[:, 0]-2451545., a[:, 1:4], a[:, 4:7], extrapolate=False)


def frozen_jpl(name):
    return state_spline(np.loadtxt(INPUTS / f"jpl_{name}.csv", delimiter=",", skiprows=1))
