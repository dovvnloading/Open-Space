"""Rank highest-visit main/bright tile centers using only catalogue quantities."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from pathlib import Path
import urllib.request
import xml.etree.ElementTree as ET

from astropy import units as u
from astropy_healpix import HEALPix
import fitsio
import numpy as np

from desi_shard import BASE, digest, obtain, product_groups

EXPOSURES = "dr1/spectro/redux/iron/exposures-iron.csv"
EXPOSURES_SHA256 = "f49d3438f0c44f4f176d8adc35bb566aa45e4615c1410001e2c78cd98ca5e2ab"
PREFIX_LIMIT_BYTES = 800_000_000


def rank(cache):
    obtain(BASE + EXPOSURES, cache / "exposures-iron.csv", EXPOSURES_SHA256)
    groups = defaultdict(list)
    with (cache / "exposures-iron.csv").open() as f:
        for row in csv.DictReader(f):
            if (row["SURVEY"], row["PROGRAM"]) == ("main", "bright"):
                groups[row["TILEID"]].append(row)
    multi = {tile: rows for tile, rows in groups.items() if len({r["NIGHT"] for r in rows}) > 1}
    max_exp = max(map(len, multi.values()))
    tiles = [(t, rows) for t, rows in multi.items() if len(rows) == max_exp]
    pix = HEALPix(nside=64, order="nested")

    def examine(pair):
        tile, visits = pair
        x = visits[0]
        ra, dec = float(x["TILERA"]), float(x["TILEDEC"])
        hp = int(pix.lonlat_to_healpix(ra * u.deg, dec * u.deg))
        redux, vac = product_groups("main", "bright", hp)
        prefix = redux[0] + "/" + redux[1][0]
        with urllib.request.urlopen(BASE + "?list-type=2&prefix=" + prefix, timeout=60) as reply:
            data = ET.fromstring(reply.read())
        size = int(data.find(".//{*}Size").text)
        name = vac[1][0]
        obtain(BASE + vac[0] + "/" + name, cache / name)
        rv = fitsio.read(str(cache / name), ext="RVTAB")
        fib = fitsio.read(str(cache / name), ext="FIBERMAP")
        ids = rv["TARGETID"]
        good = (rv["RVS_WARN"] == 0) & rv["SUCCESS"] & (rv["RR_SPECTYPE"] == "STAR") & (fib["FIBERSTATUS"] == 0)
        n = sum(np.sum(ids == t) > 1 and np.all(good[ids == t]) for t in np.unique(ids))
        return int(tile), hp, ra, dec, size, n, digest(cache / name)

    with ThreadPoolExecutor(max_workers=4) as pool:
        return sorted(pool.map(examine, tiles), key=lambda x: (-x[5], x[1]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = rank(args.cache)
    with args.output.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tileid", "healpix", "tilera", "tiledec", "spectra_bytes",
                    "rv_prefilter_repeats", "rvtab_sha256"])
        w.writerows(result)
    print("Selected:", next(r for r in result if r[4] <= PREFIX_LIMIT_BYTES))
