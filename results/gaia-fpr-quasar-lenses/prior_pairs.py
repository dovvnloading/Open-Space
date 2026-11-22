"""Supplemental vetoes for previously discussed quasar pairs and HSC lenses."""
import hashlib
import io
import json
import re
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from scipy.spatial import cKDTree


def vectors(ra, dec):
    a, d = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.column_stack((np.cos(a)*np.cos(d), np.sin(a)*np.cos(d), np.sin(d)))


def extra_catalogues(out):
    out = Path(out)
    rows, receipts = [], []

    def get(url, params=None):
        res = requests.get(url, params=params, timeout=(8, 45))
        res.raise_for_status()
        receipts.append(dict(url=res.url, size=len(res.content),
                             sha256=hashlib.sha256(res.content).hexdigest()))
        return res

    def numeric(frame, ra, dec, name=None):
        x = pd.DataFrame({"ra": pd.to_numeric(frame[ra], errors="coerce"),
                          "dec": pd.to_numeric(frame[dec], errors="coerce")})
        x["name"] = frame[name].astype(str) if name else frame.index.astype(str)
        return x[x.ra.between(0, 360) & x.dec.between(-90, 90)].copy()

    def vizier(table):
        for host in ["vizier.cds.unistra.fr", "vizier.cfa.harvard.edu"]:
            try:
                r = get("https://" + host + "/viz-bin/asu-tsv", {
                    "-source": table, "-out.max": "unlimited", "-out": "**",
                    "-out.add": "_RAJ2000,_DEJ2000", "-oc.form": "d"})
                df = pd.read_csv(io.BytesIO(r.content), sep="\t", comment="#", dtype=str)
                df.columns = df.columns.str.strip()
                for ra, dec in [("_RAJ2000", "_DEJ2000"), ("RAJ2000", "DEJ2000"),
                                ("RAdegA", "DEdegA"), ("RA_ICRS", "DE_ICRS")]:
                    if ra in df and dec in df:
                        return numeric(df, ra, dec, next((x for x in ["Name","ID","recno"] if x in df),None))
                raise ValueError("No coordinate columns in " + str(list(df)[:30]))
            except (ValueError, KeyError, requests.RequestException) as exc:
                receipts.append(dict(table=table, host=host, error=str(exc)))
        raise ValueError("No readable coordinate table")

    for tag, table, expected in [
        ("Chen et al. 2026 MGQPC", "J/A+A/707/A30/mgqpc", 4112),
        ("SuGOHI VI", "J/A+A/642/A148/sugohi", 2000),
    ]:
        try:
            records = vizier(table)
            receipts.append(dict(catalog=tag, rows=len(records), expected=expected, complete=len(records) == expected))
            records["catalog"] = tag
            rows.append(records)
        except Exception as exc:
            receipts.append(dict(catalog=tag, error=str(exc),complete=False))

    # The authors maintain a version extending the fixed 2020 CDS release.
    tag = "SuGOHI current compilation"
    try:
        page = get("https://www-utap.phys.s.u-tokyo.ac.jp/~oguri/sugohi/")
        index = BeautifulSoup(page.content, "html.parser")
        link = next(a["href"] for a in index.select("a[href]") if "download the full" in a.get_text(" ",strip=True).lower())
        data = get(urljoin(page.url, link))
        records = None
        for sep in [",", r"\s+"]:
            frame = pd.read_csv(io.BytesIO(data.content), sep=sep, engine="python",dtype=str)
            frame.columns = [re.sub("[^a-z0-9]","",str(c).lower()) for c in frame.columns]
            if {"ra","dec"}.issubset(frame.columns):
                records = numeric(frame, "ra", "dec", "name" if "name" in frame else None)
                break
        if records is None or len(records) < 2000:
            raise ValueError("Live catalogue format or population not validated")
        receipts.append(dict(catalog=tag, rows=len(records),complete=True))
        records["catalog"] = tag
        rows.append(records)
    except Exception as exc:
        receipts.append(dict(catalog=tag,error=str(exc),complete=False))

    (out / "supplemental_receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")
    result = pd.concat(rows,ignore_index=True) if rows else pd.DataFrame(columns=["ra","dec","name","catalog"])
    result.to_csv(out / "supplemental_positions.csv.gz", index=False)


def matches(components, out):
    refs = pd.read_csv(Path(out) / "supplemental_positions.csv.gz")
    tree = cKDTree(vectors(refs.ra, refs.dec))
    radius = 2 * np.sin(np.deg2rad(6 / 3600) / 2)
    result = []
    for row, nearby in zip(components.itertuples(),tree.query_ball_point(vectors(components.ra,components.dec),radius)):
        for i in nearby:
            ref = refs.iloc[i]
            result.append(dict(system_id=str(row.system_id), catalog=ref.catalog, name=ref["name"], ra=ref.ra, dec=ref.dec))
    result = pd.DataFrame(result,columns=["system_id","catalog","name","ra","dec"]).drop_duplicates()
    result.to_csv(Path(out) / "supplemental_matches.csv",index=False)
    return set(result.system_id)
