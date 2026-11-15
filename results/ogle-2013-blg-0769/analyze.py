"""Reproduce the OGLE position search and KMTNet recurrence evidence."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import html
from itertools import combinations
import json
from pathlib import Path
import re
import tarfile
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
YEARS = [*range(2011, 2020), *range(2022, 2027)]
EWS = "https://ogle.astrouw.edu.pl/ogle4/ews"
PHOT = "https://www.astrouw.edu.pl/ogle/ogle4/ews"
KMT = "https://kmtnet.kasi.re.kr/ulens/event"
POSITIONS = {
    "scene_A": "BLG03M0301.081833", "scene_B": "BLG03M0301.093760",
    "east": "BLG03M0301.032104", "west": "BLG03M0301.095181",
    "north": "BLG03M0301.035278",
}
INPUTS = {f"ews-{y}.html": f"{EWS}/{y}/ews.html" for y in YEARS}
for year, num in [(2013, "0769"), (2016, "0935")]:
    INPUTS[f"{year}-blg-{num}.dat"] = f"{PHOT}/{year}/blg-{num}/phot.dat"
    INPUTS[f"{year}-blg-{num}.html"] = f"{EWS}/{year}/blg-{num}.html"
for year, name, event, bundle in [(2016, "0228", "KB160228", "kmtp"), (2023, "0713", "KB230713", "kmt23p")]:
    INPUTS[bundle + ".gz"] = f"{KMT}/{year}/data/{event}/pysis/pysis.tar.gz"
    INPUTS[f"kmt-event-{year}.html"] = f"{KMT}/{year}/view.php?event=KMT-{year}-BLG-{name}"
for lcid in POSITIONS.values():
    INPUTS[f"{lcid}.json"] = f"https://kmtnet.kasi.re.kr/search/?action=get_lc&lcid={lcid}"
INPUTS["dn.dat"] = "https://ogle.astrouw.edu.pl/ogle/ogle4/OCVS/CV/tab1.dat"


def inputs(cache, download):
    expected = json.loads((HERE / "sources.json").read_text())
    def check(item):
        name, url = item
        path = cache / name
        if download and not path.exists():
            with urlopen(url, timeout=120) as response:
                body = response.read()
            path.write_bytes(body)
        if not path.exists():
            raise FileNotFoundError(f"Missing {name}. Supply --download to recover the input.")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if expected[name]["sha256"] != digest:
            raise ValueError(f"Archive changed for {name}. Compare with the capture before reanalysing.")
    with ThreadPoolExecutor(3) as pool:
        list(pool.map(check, INPUTS.items()))


def angle(s, hours=False):
    h, m, sec = [float(v) for v in s.split(":")]
    return np.sign(h) * (abs(h) + m/60 + sec/3600) * (15 if hours else 1)


def vectors(coords):
    ra, dec = np.deg2rad(np.array(coords).T)
    return np.stack([np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)], axis=1)


def season(path):
    events = []
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", path.read_text(errors="replace"), re.I | re.S):
        vals = [html.unescape(re.sub(r"<[^>]*>", "", v)).strip()
                for v in re.split(r"<td[^>]*>", row, flags=re.I)[1:]]
        for j, name in enumerate(vals):
            if re.fullmatch(r"[0-9]{4}-BLG-\d+", name):
                ra, dec = vals[j+3:j+5]
                events.append((name, angle(ra, True), angle(dec), vals[j+1], vals[j+2], vals[j+5]))
                break
    return events


def screen(cache, out):
    seasons = [season(cache / f"ews-{year}.html") for year in YEARS]
    coords = [vectors([(r[1], r[2]) for r in group]) for group in seasons]
    trees = [cKDTree(p) for p in coords]
    radii = np.array([0.1001, 0.3, 1, 3, 5, 10, 15, 20, 30, 60])
    chords = 2*np.sin(np.deg2rad(radii/3600)/2)
    counts = np.zeros(len(radii), dtype=int)
    closest = []
    for i, j in combinations(range(len(YEARS)), 2):
        counts += trees[i].count_neighbors(trees[j], chords)
        for k, neighbors in enumerate(trees[i].query_ball_tree(trees[j], chords[2])):
            for l in neighbors:
                sep = np.rad2deg(2*np.arcsin(np.linalg.norm(coords[i][k]-coords[j][l])/2))*3600
                a, b = seasons[i][k], seasons[j][l]
                closest.append((a[0], b[0], sep, a[3]+":"+a[4], b[3]+":"+b[4], a[5], b[5]))
    pd.DataFrame(sorted(closest, key=lambda v: v[2]), columns=["first","second","arcsec",
                 "field_star_first","field_star_second","peak_hjd_first","peak_hjd_second"]).to_csv(
                     out / "close_pairs.csv", index=False, float_format="%.6f")
    dns = [line.split() for line in (cache / "dn.dat").read_text().splitlines() if line.strip()]
    dncoords = vectors([(angle(row[1], True), angle(row[2])) for row in dns])
    target = vectors([(angle("18:05:05.39", True), angle("-27:49:21.7"))])[0]
    sep, idx = cKDTree(dncoords).query(target)
    return {"years": dict(zip(YEARS, map(len, seasons))), "total_alerts": sum(map(len,seasons)),
            "cross_year_counts": dict(zip(map(str, radii), map(int, counts))),
            "expected_unrelated_within_0.1_arcsec": {
                f"{radii[i]:g}-{radii[j]:g} arcsec": round((int(counts[j])-int(counts[i]))*0.1**2 /
                    (radii[j]**2-radii[i]**2), 5) for i,j in [(3,6),(5,7),(7,9)]},
            "OCVS_dwarf_nova_rows": len(dns), "nearest_OCVS_DN": dns[idx][0],
            "nearest_OCVS_DN_arcsec": round(float(np.rad2deg(2*np.arcsin(sep/2))*3600), 2)}


def nights(frame, center, group, unit):
    frame = frame.copy()
    frame["night"] = np.floor(frame.time).astype(int)
    bins = frame.groupby("night").value.agg(["median", "size"]).reset_index()
    bins = bins[(abs(bins.night-center)<=40) & (bins["size"]>=3)].copy()
    bins.columns = ["hjd_minus_2400000_floor", "median", "n"]
    base = bins.loc[abs(bins.iloc[:,0]-center)>=15, "median"]
    offset = base.median()
    sigma = 1.4826*(base-offset).abs().median()
    bins["baseline"] = offset
    bins["baseline_mad_scale"] = sigma
    bins["contrast"] = (bins["median"]-offset)/sigma
    bins["site_or_extraction"] = group
    bins["unit"] = unit
    bins["center"] = center
    return bins


def measurements(cache, out):
    groups = []
    for year, num, c in [(2013,"0769",56431), (2016,"0935",57536)]:
        a = np.loadtxt(cache / f"{year}-blg-{num}.dat")
        g = pd.DataFrame({"time":a[:,0]-2400000,"value":a[:,1]})
        groups.append(nights(g,c,f"OGLE_{year}","I_mag_on_artificial_baseline"))
    for filename, year, center in [("kmtp.gz",2016,57536),("kmt23p.gz",2023,60070)]:
        with tarfile.open(cache / filename) as archive:
            for site in ["A", "C", "S"]:
                a = np.loadtxt(archive.extractfile(f"KMT{site}03_I.pysis"))
                a = a[(a[:,2]>0)&(a[:,5]<9)]
                offset = 50000 if year==2016 else -2400000
                g = pd.DataFrame({"time":a[:,0]+offset,"value":-a[:,1]})
                groups.append(nights(g,center,f"pySIS_{year}_{site}","sign_reversed_difference_flux"))
    for label, lcid in POSITIONS.items():
        rows = pd.DataFrame(json.loads((cache / f"{lcid}.json").read_bytes())["rows"])
        good = rows[(rows.err>0)&(rows.secz<2)&rows.see.between(3,7)].copy()
        good = good.rename(columns={"date":"time","flux":"value"})
        for site, group in good.groupby("site"):
            for center in [57536,60070,60525]:
                groups.append(nights(group,center,label+"_"+site,"difference_flux"))
    table = pd.concat(groups,ignore_index=True)
    table.to_csv(out / "nightly.csv",index=False,float_format="%.5f")
    target = table[(table.center==60525)&(abs(table.hjd_minus_2400000_floor-60525)<=5)]
    site_results = target.groupby("site_or_extraction").contrast.agg(["count","median","max"])
    site_results["n_above_three"] = target.groupby("site_or_extraction").contrast.apply(lambda z:(z>3).sum())
    site_results.to_csv(out / "2024_checks.csv",float_format="%.3f")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    args.output.mkdir(parents=True,exist_ok=True)
    inputs(args.cache, args.download)
    summary = screen(args.cache, args.output)
    measurements(args.cache, args.output)
    (args.output / "screen.json").write_text(json.dumps(summary, indent=2)+"\n")


if __name__ == "__main__":
    main()
