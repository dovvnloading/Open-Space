"""Reproduce the OGLE cross-release screen; write the complete discordant-pair table."""
import argparse
import csv
import hashlib
import io
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np

PHASE = {3: "ogle3/OIII-CVS", 4: "ogle4/OCVS"}
URL = "https://www.astrouw.edu.pl/ogle/{}/{}/rrlyr/{}"


def fetch(url, cache):
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / hashlib.sha256(url.encode()).hexdigest()
    if not path.exists():
        for attempt in range(3):
            try:
                with urlopen(url, timeout=35) as response:
                    data = response.read()
                path.write_bytes(data)
                break
            except OSError:
                if attempt == 2:
                    raise
    return path.read_bytes()


def get(phase, cloud, filename, cache):
    return fetch(URL.format(PHASE[phase], cloud, filename), cache)


def catalog(phase, cloud, filename, cache):
    return {v[0]: v[1:] for line in get(phase, cloud, filename, cache).decode().splitlines()
            if (v := line.split())}


def photo(phase, name, band, cache):
    cloud = name.split("-")[1].lower()
    return np.loadtxt(io.BytesIO(get(phase, cloud, f"phot/{band}/{name}.dat", cache)))


def sexagesimal(s):
    h, m, sec = map(float, s.lstrip("-+").split(":"))
    return (-1 if s.startswith("-") else 1) * (h + m / 60 + sec / 3600)


def coordinates(fields):
    return np.array([15 * sexagesimal(fields[0]), sexagesimal(fields[1])])


def model(t, periods, harmonics):
    columns = [np.ones(len(t))]
    for period in periods:
        angle = 2 * np.pi * (t - np.mean(t)) / period
        for k in range(1, harmonics + 1):
            columns.extend([np.sin(k * angle), np.cos(k * angle)])
    return np.column_stack(columns)


def seasonal_strength(data, periods):
    """Rank by typical seasonal semi-amplitudes, with clipping for bad epochs."""
    if len(data) < 24:
        return None
    cut = np.floor((data[:, 0] - 1726.5) / 365.25).astype(int)
    amps = []
    for year in np.unique(cut):
        t, y, e = data[cut == year].T
        if len(t) < 24:
            continue
        X = model(t, periods, 2)
        weight = 1 / np.maximum(e, 0.015)
        keep = np.ones(len(t), dtype=bool)
        for _ in range(3):
            b = np.linalg.lstsq(X[keep] * weight[keep, None], y[keep] * weight[keep], rcond=None)[0]
            resid = y - X @ b
            scale = 1.4826 * np.median(abs(resid[keep] - np.median(resid[keep])))
            keep = abs(resid - np.median(resid)) < 4 * np.hypot(scale, e)
        amps.append(np.hypot(b[1::4], b[2::4]))
    return np.median(amps, axis=0) if amps else None


def simplify(kind):
    return {"RRe": "RRc", "aRRd": "RRd"}.get(kind, kind)


def screen(cache):
    selected, counts = [], {}
    for cloud in ["lmc", "smc"]:
        old, new = [catalog(p, cloud, "ident.dat", cache) for p in (3, 4)]
        parameters = {(p, kind): catalog(p, cloud, kind + ".dat", cache)
                      for p in (3, 4) for kind in ["RRab", "RRc", "RRd", "RRe" if p == 3 else "aRRd"]}
        shared = sorted(old.keys() & new.keys())
        inside = 0
        for name in shared:
            o, n = old[name], new[name]
            i = next(i for i, v in enumerate(o) if v.startswith("RR") or v == "aRRd")
            oldtype, newtype = o[i], n[0]
            xy = coordinates(o[i + 1:i + 3])
            xy2 = coordinates(n[1:3])
            sep = np.linalg.norm((xy - xy2) * [math.cos(math.radians(xy2[1])), 1]) * 3600
            if sep > 1:
                continue
            inside += 1
            if simplify(oldtype) == simplify(newtype):
                continue
            a, b = parameters[3, oldtype][name], parameters[4, newtype][name]
            pair = [a[2], b[2]]
            if simplify(oldtype) == "RRd":
                pair = [a[2], a[10]]
            if simplify(newtype) == "RRd":
                pair = [b[2], b[10]]
            short, long = sorted(map(float, pair))
            selected.append(dict(name=name, old=oldtype, new=newtype, sep_arcsec=round(sep, 4),
                                 short=short, long=long))
        counts[cloud] = dict(shared_id=len(shared), within_arcsec=inside)

    def measure(row):
        periods = (row["short"], row["long"])
        if not 0.70 < periods[0] / periods[1] < 0.77:
            return dict(row, status="outside radial-period ratio")
        d3, d4 = [photo(p, row["name"], "I", cache) for p in (3, 4)]
        d4 = d4[d4[:, 0] > np.max(d3[:, 0])]
        a, b = [seasonal_strength(d, periods) for d in [d3, d4]]
        if a is None or b is None:
            return dict(row, status="insufficient seasons", n3=len(d3), n4=len(d4))
        ratios = [s[0] / max(s[1], .004) for s in (a, b)]
        return dict(row, status="measured", n3=len(d3), n4=len(d4), o3=a[0], f3=a[1],
                    o4=b[0], f4=b[1], score=abs(np.log10(ratios[1] / ratios[0])))
    with ThreadPoolExecutor(8) as executor:
        measured = list(executor.map(measure, selected))
    measured.sort(key=lambda row: -row.get("score", -1))
    return counts, measured


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("raw-cache"))
    parser.add_argument("--out", type=Path, default=Path("screen.csv"))
    args = parser.parse_args()
    counts, rows = screen(args.cache)
    fields = "name old new sep_arcsec short long status n3 n4 o3 f3 o4 f4 score".split()
    with args.out.open("w", newline="") as dest:
        writer = csv.DictWriter(dest, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(counts, "classification changes", len(rows), "measured", sum(r["status"] == "measured" for r in rows))
