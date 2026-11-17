from lxml import html
from astropy.coordinates import SkyCoord
import astropy.units as u
from pathlib import Path
import pandas as pd, numpy as np, json
from urllib.request import Request, urlopen
from urllib.parse import urlencode

ref = []
r = html.parse("jades.html")
for fid in ["S3.T3", "S5.T8"]:
    for tr in r.xpath(f'//figure[@id="{fid}"]//tr')[1:]:
        c = [v.text_content().strip() for v in tr.xpath("./th|./td")]
        p = (
            [c[0:3], c[4:7]]
            if fid == "S3.T3"
            else [["JADES other variable", c[0], c[1]]]
        )
        for entry in p:
            if len(entry) != 3:
                continue
            name, ra, dec = entry
            if ":" not in ra:
                continue
            pos = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            ref.append(
                dict(
                    catalog="DeCoursey2024_" + fid,
                    designation=name,
                    ra=pos.ra.deg,
                    dec=pos.dec.deg,
                )
            )
if Path("fink.json").exists():
    for o in json.load(open("fink.json")):
        ref.append(
            dict(
                catalog="Fink_TNS",
                designation=o["f:fullname"],
                ra=float(o["f:ra"]),
                dec=float(o["f:declination"]),
            )
        )
if Path("fink_region.csv").exists():
    for r in pd.read_csv("fink_region.csv").itertuples():
        ref.append(dict(catalog="Fink_TNS", designation=r.name, ra=r.ra, dec=r.dec))
for path, cat in [
    ("cdfs.csv", "CDFS7Ms"),
    ("gaia.csv", "GaiaDR3"),
    ("smiles.csv", "Lyu2024_SMILES"),
    ("poulain.csv", "Pouliasis2019"),
    ("lyu2022.csv", "Lyu2022_AGN"),
]:
    if not Path(path).exists():
        continue
    df = pd.read_csv(path)
    for r in df.itertuples():
        ref.append(dict(catalog=cat, designation=str(r.name), ra=r.ra, dec=r.dec))
obj = pd.DataFrame(ref)
obj.to_csv("reference_positions.csv", index=False)
sites = pd.concat(
    [
        pd.read_csv("sites.csv"),
        pd.read_csv("onsets.csv")[["site", "ra", "dec"]],
        pd.read_csv("blind_pairs.csv")[["site", "ra", "dec"]],
    ]
)
rr = []
for cat, g in obj.groupby("catalog"):
    p = SkyCoord(g.ra.values, g.dec.values, unit="deg")
    q = SkyCoord(sites.ra.values, sites.dec.values, unit="deg")
    i, d, _ = q.match_to_catalog_sky(p)
    d = d.arcsec
    rr.extend(
        dict(site=s, catalog=cat, name=n, sep_arcsec=v)
        for s, n, v in zip(sites.site.values, g.designation.values[i], d)
    )
a = pd.DataFrame(rr)
a.to_csv("nearest_reference.csv", index=False)
print(
    "references",
    len(obj),
    "matches <1.5 arcsec",
    a[a.sep_arcsec < 1.5].round(3).to_string(index=False),
)
