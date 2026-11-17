from ztf_meta import c, HERE, CACHE
from query import acquire
from astropy.time import Time
from matplotlib.path import Path as Polygon
from concurrent.futures import ThreadPoolExecutor
import numpy as np, csv, json


def pick(nm):
    st = c.circular_family(nm)
    pool = []
    for year in [2026, 2025]:
        for r in csv.DictReader((CACHE / f"ztf_{nm}_{year}.raw").open()):
            if float(r["maglimit"]) < 19 or float(r["seeing"]) > 3.6:
                continue
            t = c.epoch(Time(float(r["obsjd"]), format="jd", scale="utc").isot)
            a, d = c.other_sky(st, t)
            poly = [[float(r[f"ra{k}"]), float(r[f"dec{k}"])] for k in range(1, 5)]
            good = Polygon(poly).contains_points(np.c_[a, d])
            if good.any():
                pool.append((r, good, t))
    picks = []
    n0 = len(st[0])
    hist = np.zeros(n0)
    avail = pool.copy()
    for _ in range(4):
        if not avail:
            break

        def score(v):
            r, g, t = v
            return (
                np.sum(g * (hist + 0.5))
                * (1 + 0.23 * (float(r["maglimit"]) - 20))
                * (1 + 0.04 * (int(r["filefracday"][:4]) - 2025))
            )

        r, g, t = max(avail, key=score)
        picks.append(r)
        hist += g
        avail = [v for v in avail if abs(v[2] - t) > 2.1]
    return picks


def url(row, suffix):
    f = row["filefracday"]
    fields = [
        f,
        row["field"].zfill(6),
        row["filtercode"],
        "c" + row["ccdid"].zfill(2),
        row["imgtypecode"],
        "q" + row["qid"],
    ]
    return (
        "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/"
        + f[:4]
        + "/"
        + f[4:8]
        + "/"
        + f[8:]
        + "/ztf_"
        + "_".join(fields)
        + "_"
        + suffix
    )


if __name__ == "__main__":
    choose = {nm: pick(nm) for nm in ["Chen1", "Chen2", "RR"]}
    (HERE / "ztf_selection.json").write_text(json.dumps(choose, indent=2) + "\n")
    items = {}
    for nm, rows in choose.items():
        print(nm, [(x["filefracday"], x["pid"]) for x in rows], flush=True)
        for r in rows:
            for name, suffix in [
                ("diff", "scimrefdiffimg.fits.fz"),
                ("psf", "diffimgpsf.fits"),
                ("mask", "mskimg.fits"),
            ]:
                items["ztf_" + r["pid"] + "_" + name] = url(r, suffix)
    with ThreadPoolExecutor(8) as e:
        list(e.map(lambda it: acquire(*it), items.items()))
