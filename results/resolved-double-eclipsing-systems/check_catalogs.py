"""Recompute nearest literature positions for the auditable 15-pair cohort."""
import argparse
import csv
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.spatial import cKDTree

CATALOGS = {
    "hartman": "J/ApJ/934/72",
    "msc": "J/ApJS/235/6",
    "vaessen": "J/A+A/682/A164",
    "kostov": "J/AJ/171/29",
    "kostov1": "J/ApJS/259/66",
    "kostov2real": "J/MNRAS/527/3995",
}


def vectors(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)],axis=-1)


def read_points(content):
    header = []
    coords = []
    for line in content.splitlines():
        cells = line.split("\t")
        if line.startswith("#Table"):
            header = []
        if cells[0] == "recno":
            header = cells
        if cells[0].strip().isdigit():
            row = dict(zip(header,cells))
            if "RAJ2000" in row:
                coords.append((float(row["RAJ2000"]),float(row["DEJ2000"])))
            elif "RA_ICRS" in row:
                coords.extend([(float(row["RA_ICRS"]),float(row["DE_ICRS"])),(float(row["RA2deg"]),float(row["DE2deg"]))])
            elif "_RA" in row:
                coords.append((float(row["_RA"]),float(row["_DE"])))
    if not coords:
        raise ValueError("No catalog coordinates")
    return coords


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("pairs",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    rows = list(csv.DictReader(args.pairs.open()))
    refs = set(Path(__file__).with_name("fezenko_ids.txt").read_text().split())
    results = {r["source_id1"]:{"fezenko2022_id": bool(refs & {r["source_id1"],r["source_id2"]})} for r in rows}
    for key,cat in CATALOGS.items():
        params = urlencode({"-source":cat,"-out":"**","-out.max":50000})
        with urlopen("https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"+params,timeout=120) as f:
            positions = read_points(f.read().decode())
        tree = cKDTree(vectors(*np.array(positions).T))
        for r in rows:
            stars = [vectors(float(r["ra"+i]),float(r["dec"+i])) for i in "12"]
            delta = tree.query(stars)[0].min()
            results[r["source_id1"]][key+"_min_arcsec"]=round(np.degrees(2*np.arcsin(delta/2))*3600,1)
    args.output.write_text(json.dumps(results,indent=2)+"\n")
