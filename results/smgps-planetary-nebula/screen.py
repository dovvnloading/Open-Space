"""Cross-match published SMGPS candidate centers against MASH-I and MASH-II."""

import csv
import json
import urllib.request
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
BASE = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
URLS = {
    "snr.tsv": BASE
    + "-source=J/A%2bA/693/A247/meerkatc&-out=Name,GLON,GLAT,Rad,Rel&-out.max=unlimited",
    "mash1.tsv": BASE
    + "-source=V/127A/mash1&-out=PNG,n_PNG,Name,GLon,GLat,MajDiam,MinDiam&-out.max=unlimited",
    "mash2.tsv": BASE
    + "-source=V/127A/mash2&-out=PNG,n_PNG,Name,GLon,GLat,MajDiam,MinDiam&-out.max=unlimited",
}


def get(name, url):
    CACHE.mkdir(exist_ok=True)
    p = CACHE / name
    if not p.exists():
        with urllib.request.urlopen(url, timeout=120) as r:
            data = r.read()
        p.write_bytes(data)
    return p


def rows(name):
    lines = get(name, URLS[name]).read_text().splitlines()
    rows = [line for line in lines if line.strip() and not line.startswith("#")]
    return list(csv.DictReader([rows[0]] + rows[3:], delimiter="\t"))


def compare():
    # SMGPS repeats one row verbatim. Spatial eligibility requires a measured
    # size on both surveys and a separation within one fifth of the radio radius.
    objects = {s["Name"].strip(): s for s in rows("snr.tsv")}
    pns = rows("mash1.tsv") + rows("mash2.tsv")
    pns = [
        p
        for p in pns
        if all(p[k].strip() for k in ["GLon", "GLat", "MajDiam", "MinDiam"])
    ]
    stars = SkyCoord(
        [float(p["GLon"]) for p in pns] * u.deg,
        [float(p["GLat"]) for p in pns] * u.deg,
        frame="galactic",
    )
    targets = list(objects.values())
    glon = np.array([float(s["GLON"]) for s in targets])
    glat = np.array([float(s["GLAT"]) for s in targets])
    radius = np.array([float(s["Rad"]) * 60 for s in targets])
    diam = np.array([(float(p["MajDiam"]) * float(p["MinDiam"])) ** 0.5 for p in pns])

    def hits(shift):
        from astropy.coordinates import search_around_sky

        cc = SkyCoord((glon + shift) * u.deg, glat * u.deg, frame="galactic")
        ii, jj, dist, _ = search_around_sky(cc, stars, np.max(radius) * 0.2 * u.arcsec)
        ans = []
        for i, j, d in zip(ii, jj, dist.arcsec):
            scale = diam[j] / (2 * radius[i])
            if d < 0.2 * radius[i] and 0.6 <= scale <= 1.5:
                ans.append(
                    dict(
                        radio=targets[i]["Name"].strip(),
                        pn=pns[j]["Name"].strip(),
                        classification=pns[j]["n_PNG"],
                        offset_arcsec=round(d, 4),
                        size_ratio=round(scale, 4),
                    )
                )
        return ans

    real = hits(0)
    # Translations of the candidate longitudes preserve the narrow latitude distribution.
    shifts = [i / 10 for i in range(-50, 51) if abs(i) >= 3]
    summary = dict(
        radio_rows=len(rows("snr.tsv")),
        distinct_radio=len(targets),
        mash1_rows=len(rows("mash1.tsv")),
        mash2_rows=len(rows("mash2.tsv")),
        eligible_optical=len(pns),
        match=real,
        criteria="offset < 0.2 radio radius; equivalent PN diameter / radio diameter 0.6 to 1.5 inclusive",
        shifted_counts={str(x): len(hits(x)) for x in shifts},
    )
    (HERE / "catalog-test.json").write_text(json.dumps(summary, indent=2) + "\n")
    print({k: v for k, v in summary.items() if k != "shifted_counts"})
    print("Shifted total", sum(summary["shifted_counts"].values()))


if __name__ == "__main__":
    compare()
