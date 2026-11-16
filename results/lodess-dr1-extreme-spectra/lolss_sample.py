"""Screen compact LoLSS DR1 sources and measure signed LoDeSS flux at their positions."""
import argparse
import hashlib
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.spatial import cKDTree

from comparison.match_catalogues import unit_xyz
from release.photometry import RadioImage

FREQUENCIES = np.arange(44, 65, 4)
BASE_URL = "https://lofar-surveys.org/public/lolss/catalogues/LoLSS_DR1_v1.1"


def read_table(path):
    table = Table.read(path).to_pandas()
    for col in table.select_dtypes("object"):
        table[col] = table[col].str.decode("ascii").str.strip()
    return table


def nearest(query, reference):
    tree = cKDTree(unit_xyz(reference.RA, reference.DEC))
    dist, ix = tree.query(unit_xyz(query.RA, query.DEC), k=2)
    return np.degrees(2 * np.arcsin(dist / 2)) * 3600, ix


def fit_slope(row):
    flux = row[[f"s{f}" for f in FREQUENCIES]].to_numpy(dtype=float)
    err = row[[f"e{f}" for f in FREQUENCIES]].to_numpy(dtype=float)
    ok = np.isfinite(flux + err) & (flux > 0) & (err > 0)
    if ok.sum() < 5:
        return np.full(3, np.nan)
    # A common flux scale cancels for the slope. Adopt a 4% differential floor.
    w = 1 / ((err[ok] / flux[ok])**2 + .04**2)
    design = np.column_stack((np.ones(ok.sum()), np.log(FREQUENCIES[ok] / 54)))
    cov = np.linalg.inv((design.T * w) @ design)
    par = cov @ (design.T @ (w * np.log(flux[ok])))
    chi = sum(w * (np.log(flux[ok]) - design @ par)**2)
    error = np.sqrt(cov[1, 1] * max(1, chi / (ok.sum() - 2)))
    return par[1], error, chi


def run(release, cache, out):
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    receipts = []
    for name, suffix in [("lolss", "")] + [(f"chan{i}", f"-chan{i:02}") for i in range(6)]:
        url = BASE_URL + suffix + ".srl.fits"
        path = cache / (name + ".fits")
        if not path.exists():
            tmp = path.with_suffix(".tmp")
            with urlopen(url, timeout=120) as response:
                tmp.write_bytes(response.read())
            Table.read(tmp, format="fits")
            tmp.replace(path)
        receipts.append((url, path.stat().st_size, hashlib.sha256(path.read_bytes()).hexdigest()))
    pd.DataFrame(receipts, columns=["url", "bytes", "sha256"]).to_csv(out / "lolss_inputs.csv", index=False)
    src = read_table(cache / "lolss.fits")
    f = pd.read_csv(release / "fields.csv").query("quality == 'good'").rename(
        columns={"ra_deg": "RA", "dec_deg": "DEC"})
    d, ix = nearest(src, f)
    src["field"] = f.field.to_numpy()[ix[:, 0]]
    src["radius"] = d[:, 0] / 3600
    d, _ = nearest(src, src)
    src["nn"] = d[:, 1]
    src = src.loc[(src.radius <= 3.5) & (src.Total_flux >= 100) & (src.nn >= 60)
                  & (src.Total_flux / src.Peak_flux < 1.35) & (src.S_Code == "S")].copy()
    lodess = read_table(release / "sources.fits")
    d, ix = nearest(src, lodess)
    src["lo_idx"] = np.where(d[:, 0] < 30, ix[:, 0], -1)
    for i, freq in enumerate(FREQUENCIES):
        ch = read_table(cache / f"chan{i}.fits")
        d, ix = nearest(src, ch)
        for col, key in [("Total_flux", "s"), ("E_Total_flux", "e")]:
            src[f"{key}{freq}"] = np.where(d[:, 0] < 12, ch[col].to_numpy()[ix[:, 0]], np.nan)
    src[["slope", "eslope", "chisq"]] = np.array([fit_slope(r) for _, r in src.iterrows()])
    phot = {}
    for name, group in src.groupby("field"):
        with RadioImage(release / "images" / (name + "-I.fits")) as image:
            for i, r in group.iterrows():
                phot[i] = image.measure(r.RA, r.DEC) or {}
    src = src.join(pd.DataFrame.from_dict(phot, orient="index"))
    src.index.name = "lolrow"
    src.to_csv(out / "lolss_parent.csv", float_format="%.9g")
    Table([src.RA.to_numpy(), src.DEC.to_numpy()], names=["RA", "DEC"]).write(
        cache / "parent_positions.fits", overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ["release_dir", "cache", "output"]:
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    run(args.release_dir, args.cache, args.output)
