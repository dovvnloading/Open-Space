"""Replay the exploratory linear track searches from both public source lists."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
START = 1990
END = 2120
MAS_PER_DEG = 3_600_000


def acquire(row, cache):
    dest = cache / row["file"]
    if not dest.is_file():
        temp = dest.with_suffix(dest.suffix + ".part")
        with urlopen(row["url"], timeout=600) as response, temp.open("wb") as out:
            while buf := response.read(1024*1024):
                out.write(buf)
        temp.replace(dest)
    with dest.open("rb") as stream:
        checksum = hashlib.sha256()
        for buf in iter(lambda: stream.read(1024*1024), b""):
            checksum.update(buf)
    if checksum.hexdigest() != row["sha256"]:
        # TAP does not guarantee row order. Comparing all retained pairs and
        # the original counts checks a new export. Binary releases are fixed.
        if dest.suffix != ".csv":
            raise ValueError(f"Input release hash differs: {dest}")
        print(f"Changed TAP response bytes: {dest.name}; checking retained pairs")
    return dest


def unit_vectors(ra, dec):
    a, d = np.radians([ra, dec])
    return np.column_stack((np.cos(a)*np.cos(d), np.sin(a)*np.cos(d), np.sin(d)))


def quasars(path):
    if path.suffix == ".zip":
        rows = []
        with ZipFile(path) as archive, archive.open("milliquas.txt") as catalog:
            for line in catalog:
                if line[51:52] in (b"Q", b"A"):
                    rows.append((line[25:50].decode().strip(), float(line[:11]),
                                 float(line[12:23]), float(line[76:82].strip() or 0)))
        return pd.DataFrame(rows, columns=["id", "ra", "dec", "z"])
    with fits.open(path) as fits_file:
        data = fits_file[1].data
        usable = (data["COADD_FIBERSTATUS"] == 0) & (data["ZWARN"] == 0) & (data["Z"] > 0.1)
        columns = {"id":"TARGETID", "ra":"TARGET_RA", "dec":"TARGET_DEC", "z":"Z"}
        qso = pd.DataFrame({key:data[name][usable].astype(data[name].dtype.newbyteorder("="))
                            for key, name in columns.items()})
        return qso.drop_duplicates("id").reset_index(drop=True)


def screen(stars, sources):
    tree = cKDTree(unit_vectors(sources.ra, sources.dec))
    mid = (START + END)/2 - 2016
    half = (END - START)/2
    cosdec = np.cos(np.radians(stars.dec))
    speed = np.hypot(stars.pmra, stars.pmdec)
    midpoint = unit_vectors(stars.ra + mid*stars.pmra/(cosdec*MAS_PER_DEG),
                            stars.dec + mid*stars.pmdec/MAS_PER_DEG)
    groups = tree.query_ball_point(midpoint, np.radians((half*speed + 4000)/MAS_PER_DEG))
    left = np.repeat(np.arange(len(stars)), [len(x) for x in groups])
    right = np.concatenate([x for x in groups if len(x)]).astype(int)
    lens = stars.iloc[left].reset_index(drop=True)
    qso = sources.iloc[right].reset_index(drop=True)
    x = ((qso.ra-lens.ra+180) % 360-180)*np.cos(np.radians(lens.dec))*MAS_PER_DEG
    y = (qso.dec-lens.dec)*MAS_PER_DEG
    dt = (x*lens.pmra + y*lens.pmdec)/(lens.pmra**2 + lens.pmdec**2)
    out = lens.join(qso, rsuffix="_q")
    out["miss"] = np.hypot(x - lens.pmra*dt, y - lens.pmdec*dt)
    out["jyear"] = 2016 + dt
    out["sep16"] = np.hypot(x, y)/1000
    return out[out.jyear.between(START, END) & (out.miss < 4000)].sort_values("miss")


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("cache", type=Path, help="download directory, about 1 GB")
    cache = args.parse_args().cache
    cache.mkdir(parents=True, exist_ok=True)
    manifest = {r["file"]:r for r in json.loads((HERE/"sources.json").read_text())}
    frames = [pd.read_csv(acquire(manifest[n], cache)).dropna()
              for n in ("gaia.csv", "gcns.csv", "gcsext.csv")]
    stars = pd.concat(frames, ignore_index=True).drop_duplicates("source_id").reset_index(drop=True)
    if len(stars) != 1276314:
        raise ValueError(f"Foreground coverage changed to {len(stars)}")
    for name, filename, count in [("mq","milliquas.zip",907144), ("desi","qso.fits",1589836)]:
        sources = quasars(acquire(manifest[filename], cache))
        if len(sources) != count:
            raise ValueError(f"{name} coverage changed")
        result = screen(stars, sources)
        expected = pd.read_csv(HERE / f"screen-{name}.csv")
        for t in (result, expected):
            t["id"] = t.id.astype(str)
        keys = ["source_id", "id"]
        a = result.sort_values(keys).reset_index(drop=True)
        b = expected.sort_values(keys).reset_index(drop=True)
        pd.testing.assert_frame_equal(a[keys], b[keys])
        np.testing.assert_allclose(a[["miss","jyear"]], b[["miss","jyear"]], atol=1e-5, rtol=0)
        print(f"{name}: {len(stars)} foreground sources, {len(sources)} entries, {len(result)} retained pairs verified")


if __name__ == "__main__":
    main()
