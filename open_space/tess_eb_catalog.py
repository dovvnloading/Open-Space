"""Search the archived TESS eclipsing-binary catalog for threshold-crossing detections."""

import argparse
import csv
from hashlib import sha256
import io
import json
import math
from pathlib import Path

from .tess_archive import atomic_bytes, download


CATALOG_URL = ("https://archive.stsci.edu/hlsps/tess-ebs/"
               "hlsp_tess-ebs_tess_lcf-ffi_s0001-s0026_tess_v1.0_cat.csv")
CATALOG_HASH = "97e81f6c01716a431379ab6bd37d1182b611a485107110e0eba830f817ba7a8f"


def separation(ra, dec, ra2, dec2):
    lat1, lat2 = math.radians(dec), math.radians(dec2)
    a = (math.sin((lat2-lat1)/2)**2 + math.cos(lat1)*math.cos(lat2)*
         math.sin(math.radians(ra2-ra)/2)**2)
    return math.degrees(2*math.asin(math.sqrt(min(1.0,a))))*3600


def load_catalog(cache_dir):
    catalog = Path(cache_dir) / "tess-ebs-v1.0.csv"
    raw = catalog.read_bytes() if catalog.exists() else download(CATALOG_URL, 8_000_000)
    if sha256(raw).hexdigest() != CATALOG_HASH:
        raise ValueError("The MAST EB catalog differs from the frozen v1 snapshot")
    if not catalog.exists():
        atomic_bytes(catalog, raw)
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    return rows


def lookup(catalog_rows, tic, ra_deg, dec_deg, radius_arcsec=90):
    matches = []
    for row in catalog_rows:
        same_id = int(row["tess_id"]) == tic
        distance = separation(ra_deg, dec_deg,
                              float(row["ra"]), float(row["dec"]))
        if same_id or distance <= radius_arcsec:
            matches.append({"catalog_tic": int(row["tess_id"]), "same_id": same_id,
                            "separation_arcsec": distance, "catalog_period_days": float(row["period"])})
    return sorted(matches, key=lambda r: r["separation_arcsec"])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence", type=Path, nargs="+")
    p.add_argument("--cache", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    rows = load_catalog(args.cache)
    cases = []
    for source in args.evidence:
        evidence = json.loads(source.read_bytes())
        for target in evidence["targets"]:
            if target["role"] == "known_control" or target["disposition"] not in (
                    "two_sector_periodic_dip", "periodic_dip_with_veto"):
                continue
            cases.append({"tic_id":target["tic_id"],
                          "matched":lookup(rows,target["tic_id"],target["ra_deg"],target["dec_deg"])})
    result = {"catalog_url":CATALOG_URL, "catalog_sha256":CATALOG_HASH, "total_rows":len(rows),
              "radius_arcsec":90, "epoch_propagation":"none", "cases":cases}
    atomic_bytes(args.output,(json.dumps(result,indent=2,allow_nan=False)+"\n").encode())


if __name__=="__main__":
    main()
