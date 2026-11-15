"""Find Galactic recombination-line systems in two public FASHI H I releases."""
import argparse
import hashlib
import itertools
import json
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

C = 299792.458
# J2000 ICRS -> Galactic rotation, and standard kinematic LSR solar motion.
ROT = np.array([[-.0548755604, -.8734370902, -.4838350155],
                [.4941094279, -.4448296300, .7469822445],
                [-.8676661490, -.1980763734, .4559837762]])
SOLAR = np.array([10.27, 15.32, 7.74])
ARCSEC = 180 * 3600 / np.pi
FILES = {
    "fashi.csv": ("https://zcp521.github.io/img/Table2-FASHI_extragalactic_HI_source_catalog.csv",
                  "e0346a4bedba17630d52f364676417a0fb270777ec3bcf62d57b2efddc6712a6"),
    "dr2.csv": ("https://nadc.china-vo.org/res/file_upload/download?id=1538184",
                "9c612acac23776b99630f72330f818ca01b6e76517d0107a7872b056b703db4a"),
    "wise.csv": ("https://astro.phys.wvu.edu/wise/catalog/WISE_HII_V3.0.csv",
                 "a1c78ce39afa4e6b8f3e0ed7361d8487f0d235257a15abe5fcda6f0d267856eb"),
}
# Hydrogen Rydberg frequency, including the electron/proton reduced mass.
def rest(n, step):
    return 10973731.568157 * 299792458 / (1 + 1 / 1836.152673426) * (
        n**-2 - (n + step)**-2) / 1e6


def unit(ra, dec):
    r, d = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])


def inputs(directory):
    directory.mkdir(parents=True, exist_ok=True)
    for name, (url, sha) in FILES.items():
        path = directory / name
        if not path.exists():
            with tempfile.NamedTemporaryFile(dir=directory) as temp:
                urllib.request.urlretrieve(url, temp.name)
                if hashlib.sha256(Path(temp.name).read_bytes()).hexdigest() != sha:
                    raise ValueError(f"Release content changed: {url}")
                path.write_bytes(Path(temp.name).read_bytes())
        if hashlib.sha256(path.read_bytes()).hexdigest() != sha:
            raise ValueError(f"Unexpected checksum for {path}")


def load(directory):
    dr1 = pd.read_csv(directory / "fashi.csv", skiprows=[1]).rename(columns={
        "Name": "name", "ID_FASHI": "id", "SNR": "snr", "log10Mass": "mass"})
    dr2 = pd.read_csv(directory / "dr2.csv").rename(columns={"fashi_id": "id"})
    wise = pd.read_csv(directory / "wise.csv", encoding="latin1", low_memory=False)
    regions = wise.dropna(subset=["VLSR", "RA (J2000)", "Dec (J2000)"]).reset_index(drop=True)
    catalogs = []
    for release, catalog in enumerate([dr1, dr2], 1):
        catalog["release"] = release
        catalog = catalog.loc[catalog.W_50.between(0, 100, inclusive="neither")].copy()
        vectors = unit(catalog.ra, catalog.dec)
        catalog["lsr_correction"] = vectors @ ROT.T @ SOLAR
        catalogs.append((catalog.reset_index(drop=True), cKDTree(vectors)))
    return catalogs, regions


def matches(catalog_item, regions, offset=0., extension=False):
    catalog, tree = catalog_item
    sky = unit(regions["RA (J2000)"], regions["Dec (J2000)"])
    # The counterfactual shifts each region in Galactic longitude at fixed latitude.
    if offset:
        gal = sky @ ROT.T
        lon = np.rad2deg(np.arctan2(gal[:, 1], gal[:, 0])) + offset
        lat = np.rad2deg(np.arcsin(gal[:, 2]))
        sky = unit(lon, lat) @ ROT
    max_radius = regions.Radius.max() if extension else 180.
    pairs = tree.sparse_distance_matrix(cKDTree(sky), 2*np.sin(max_radius/(2*ARCSEC)),
                                        output_type="ndarray")
    i, j = pairs["i"], pairs["j"]
    sep = 2 * ARCSEC * np.arcsin(pairs["v"]/2)
    keep = sep <= regions.Radius.values[j]
    i, j, sep = i[keep], j[keep], sep[keep]
    result = []
    transitions = [(n, 1) for n in range(167, 171)]
    if extension:
        transitions += [(n, 2) for n in range(208, 216)]
    for n, step in transitions:
        freq = rest(n, step)
        vhel = C * (1 - catalog.freq.values[i]/freq)
        vlsr = vhel + catalog.lsr_correction.values[i]
        keep = abs(vlsr-regions.VLSR.values[j]) <= 16
        cols = ["release", "id", "name", "ra", "dec", "freq", "freq_err", "W_50",
                "snr", "distance", "mass"]
        rows = catalog.iloc[i[keep]][cols].copy()
        rows["region"] = regions.WISE_Name.values[j[keep]]
        rows["n"] = n
        rows["series"] = "alpha" if step == 1 else "beta"
        rows["rest_mhz"] = freq
        rows["sep_arcsec"] = sep[keep]
        rows["vhel_kms"] = vhel[keep]
        rows["vlsr_kms"] = vlsr[keep]
        rows["reference_kms"] = regions.VLSR.values[j[keep]]
        result.append(rows)
    return pd.concat(result, ignore_index=True)


def seeds(rows):
    selected = set()
    for region, group in rows.groupby("region"):
        if group.n.nunique() < 2:
            continue
        for a, b in itertools.combinations(group.itertuples(), 2):
            if a.n == b.n or abs(a.vlsr_kms - b.vlsr_kms) > 12:
                continue
            pair = unit([a.ra, b.ra], [a.dec, b.dec])
            if np.linalg.norm(pair[0]-pair[1])*ARCSEC <= 180:
                selected.add(region)
                break
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    inputs(args.cache)
    catalogs, regions = load(args.cache)
    initial = [seeds(matches(cat, regions)) for cat in catalogs]
    selected = regions[regions.WISE_Name.isin(set.union(*initial))].reset_index(drop=True)
    result = pd.concat([matches(cat, selected, extension=True) for cat in catalogs])
    result = result.sort_values(["release", "region", "series", "n", "id"])
    args.output.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output/"identifications.csv", index=False, float_format="%.6f")
    columns = ["WISE_Name", "All_Names", "RA (J2000)", "Dec (J2000)", "Radius", "VLSR",
               "e_VLSR", "FWHM", "Lines", "Author"]
    selected[columns].to_csv(args.output/"regions.csv", index=False)
    controls = []
    for offset in np.r_[-np.arange(.5, 10.1, .5), np.arange(.5, 10.1, .5)]:
        found = [seeds(matches(cat, regions, offset)) for cat in catalogs]
        controls.append({"delta_l_deg": offset, "dr1": len(found[0]), "dr2": len(found[1])})
    summary = {"raw_radio_rows": [41741, 156411], "width_selected_radio_rows": [
        len(cat) for cat, _ in catalogs], "regions_with_velocities": len(regions),
        "seed_regions_by_release": [sorted(x) for x in initial], "offset_controls": controls}
    (args.output/"audit.json").write_text(json.dumps(summary, indent=2)+"\n")


if __name__ == "__main__":
    main()
