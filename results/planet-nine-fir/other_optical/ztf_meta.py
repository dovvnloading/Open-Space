from query import acquire, ThreadPoolExecutor, urlencode, HERE, CACHE
from astropy.time import Time
from concurrent.futures import ThreadPoolExecutor
import sys, json, numpy as np

for path in [HERE.parent.parent, HERE.parent]:
    if (path / "other" / "candidates.py").exists():
        sys.path.insert(0, str(path / "other"))
        break
import candidates as c

RADII = {
    "Chen1": (42.676, -15.021, 1.9),
    "Chen2": (45.297, -16.711, 1.9),
    "RR": (0.0, 0.0, 0.0),
}


def run(key):
    a, d, sz = RADII[key]
    t = c.epoch("2026-09-01")
    st = c.circular_family(key)
    x, y = c.other_sky(st, t)
    if key == "RR":
        # A bounded local subset of the long ring: point midway among the first 7 headings.
        a, d, sz = float(np.median(x[:7])), float(np.median(y[:7])), 1.5
    for y in [2026, 2025]:
        where = f"obsjd BETWEEN {Time(str(y)+'-08-01').jd} AND {Time(str(y)+'-11-17').jd} AND filtercode='zr'"
        url = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?" + urlencode(
            {
                "POS": f"{a:.5f},{d:.5f}",
                "SIZE": str(sz),
                "WHERE": where,
                "RESPONSEFORMAT": "CSV",
                "columns": "filefracday,field,filtercode,ccdid,imgtypecode,qid,obsjd,maglimit,seeing,ra,dec,ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4,pid",
            }
        )
        acquire(f"ztf_{key}_{y}", url)


if __name__ == "__main__":
    with ThreadPoolExecutor(3) as e:
        list(e.map(run, RADII))
