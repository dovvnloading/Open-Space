"""Query the two observing teams' public APT project catalogs for prior transient flags."""

from pathlib import Path
from io import StringIO, BytesIO
from zipfile import ZipFile
from xml.etree.ElementTree import fromstring
from urllib.request import urlopen
from hashlib import sha256
from scipy.spatial import cKDTree
import pandas as pd, numpy as np, sys

root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
sites = pd.concat(
    [
        pd.read_csv("sites.csv"),
        pd.read_csv("onsets.csv")[["site", "ra", "dec"]],
        pd.read_csv("blind_pairs.csv")[["site", "ra", "dec"]],
    ]
)
records = []
receipts = []
for pid in [8060, 12577]:
    url = f"https://www.stsci.edu/jwst-program-info/download/jwst/apt/{pid}/"
    p = root / f"{pid}.aptx"
    data = p.read_bytes() if p.exists() else urlopen(url, timeout=60).read()
    receipts.append(dict(pid=pid, url=url, sha256=sha256(data).hexdigest()))
    t = fromstring(ZipFile(BytesIO(data)).read(f"{pid}.xml"))
    for elt in t.iter():
        if elt.tag.split("}")[-1] != "Catalog":
            continue
        children = {a.tag.split("}")[-1]: a.text for a in elt}
        txt = children.get("CatalogAsCsv")
        if not txt or "[LABEL] - Class" not in txt.splitlines()[0]:
            continue
        d = pd.read_csv(StringIO(txt))
        d = d[d["[LABEL] - Class"].fillna("").str.contains("trans", case=False)]
        if not len(d):
            continue
        xy = lambda ra, dec: np.c_[ra * np.cos(np.deg2rad(-27.79)), dec] * 3600
        dist, j = cKDTree(xy(d.RA, d.DEC)).query(xy(sites.ra, sites.dec))
        for s, dt, idx in zip(sites.site, dist, j):
            row = d.iloc[idx]
            records.append(
                dict(
                    site=s,
                    pid=pid,
                    catalog=children.get("CatalogName") or children.get("Name"),
                    n_transient_or_host=len(d),
                    target_id=row["#ID"],
                    target_label=row["[LABEL] - Class"],
                    ra=row.RA,
                    dec=row.DEC,
                    sep_arcsec=dt,
                )
            )
pd.DataFrame(receipts).to_csv("apt_inputs.csv", index=False)
r = pd.DataFrame(records)
r.to_csv("apt_nearest.csv", index=False)
print(
    r.sort_values("sep_arcsec")
    .drop_duplicates("site")
    .query("site in [93,99,185,525,829,1004,2010,2027,2028]")[
        ["site", "target_id", "target_label", "sep_arcsec"]
    ]
    .round(3)
    .to_string(index=False)
)
