#!/usr/bin/env python3
"""Freeze and verify DJA v4.4 source bytes; no silent release substitution.

Requires Python 3.10+. Network only for 'download' and 'mast-index'.
Reproduction:
  python retrieve.py download --catalog-only --cache /data/jwst-dja
  python retrieve.py download --all --cache /data/jwst-dja
  python retrieve.py freeze --cache /data/jwst-dja
  python retrieve.py mast-index --cache /data/jwst-dja
"""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import urllib.request

HERE = Path(__file__).resolve().parent
SOURCES = json.loads((HERE / "frozen_sources.json").read_text())
COLUMNS = (
    "index", "root", "file", "srcid", "objid", "ra", "dec",
    "grating", "filter", "grade", "z_best", "z_prism", "z_grating",
    "dataset", "msamet", "nfiles", "phot_dr"
)
SHARDS = {row["band"]: row for row in SOURCES["dja"]["spectral_shards"]}


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def digests(path):
    sha, md5 = hashlib.sha256(), hashlib.md5()
    count = 0
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(2**20), b""):
            sha.update(block)
            md5.update(block)
            count += len(block)
    return {"sha256": sha.hexdigest(), "md5": md5.hexdigest(), "bytes": count}


def safe_write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temp.replace(path)


def download_one(item, cache):
    name = item["filename"]
    path = cache / name
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        tmp = path.with_suffix(path.suffix + ".partial")
        try:
            req = urllib.request.Request(item["url"], headers={"User-Agent": "Open-Space-repro/1.0"})
            with urllib.request.urlopen(req, timeout=120) as inp, open(tmp, "wb") as out:
                while True:
                    data = inp.read(2**20)
                    if not data:
                        break
                    out.write(data)
            tmp.replace(path)
        except Exception:
            tmp.unlink(missing_ok=True)
            raise
    receipt = digests(path)
    receipt.update({"file": name, "url": item["url"], "verified_at_utc": utc_now()})
    if receipt["md5"].lower() != item["publisher_md5"].lower():
        raise ValueError(f"Publisher MD5 mismatch for {name}; refuse to freeze.")
    receipt["publisher_md5_matches"] = True
    return receipt


def items(all_files):
    dja = SOURCES["dja"]
    return [dja["catalog"], dja["columns"]] + (dja["spectral_shards"] if all_files else [])


def cmd_download(args):
    cache = args.cache
    ledger = []
    for item in items(args.all):
        receipt = download_one(item, cache)
        print(f"Verified {receipt['file']}: sha256={receipt['sha256']}")
        ledger.append(receipt)
        safe_write_json(cache / "download_receipts.json", {"source_record": dja_record(), "files": ledger})


def dja_record():
    return {"record_url": SOURCES["dja"]["record_url"], "version": SOURCES["dja"]["version"]}


def band_key(grating, filt):
    g, f = grating.lower().strip(), filt.lower().strip()
    return "prism" if g == "prism" else f"{g}-{f}"


def cmd_freeze(args):
    cache = args.cache
    catalog = cache / SOURCES["dja"]["catalog"]["filename"]
    if not catalog.exists():
        raise FileNotFoundError("Download verified DJA summary CSV first; cannot freeze unobserved rows.")
    if digests(catalog)["md5"] != SOURCES["dja"]["catalog"]["publisher_md5"]:
        raise ValueError("Input catalog is not the frozen v4.4 snapshot.")
    out = cache / "frozen_row_manifest.csv"
    counts = {key: 0 for key in SHARDS}
    written = 0
    with gzip.open(catalog, "rt", encoding="utf-8", newline="") as inp, \
         open(out.with_suffix(".partial"), "w", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(inp)
        missing = set(COLUMNS[1:]) - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Unexpected DJA schema: " + ", ".join(sorted(missing)))
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        for i, row in enumerate(reader):
            k = band_key(row["grating"], row["filter"])
            if k not in SHARDS:
                raise ValueError(f"Unknown grating/filter at row {i}: {k}")
            record = {col: row.get(col, "") for col in COLUMNS}
            record["index"] = i
            writer.writerow(record)
            counts[k] += 1
            written += 1
    if written != SOURCES["dja"]["expected_entries"]:
        out.with_suffix(".partial").unlink(missing_ok=True)
        raise AssertionError(f"Expected 80367 spectrum records, got {written}")
    for key, count in counts.items():
        if count != SHARDS[key]["expected_spectra"]:
            out.with_suffix(".partial").unlink(missing_ok=True)
            raise AssertionError(f"{key} count {count} != {SHARDS[key]['expected_spectra']}")
    out.with_suffix(".partial").replace(out)
    rec = digests(out)
    safe_write_json(cache / "freeze_receipt.json", {
        "generated_at_utc": utc_now(),
        "catalog_sha256": digests(catalog)["sha256"],
        "manifest_sha256": rec["sha256"],
        "rows": written,
        "counts": counts,
        "source_release": dja_record(),
    })
    print(json.dumps({"status": "FROZEN", "rows": written, "manifest_sha256": rec["sha256"]}))


def cmd_mast_index(args):
    # A *separate* public-product inventory: MAST program metadata are not extracted 1D spectra.
    try:
        from astroquery.mast import Observations
    except ImportError as exc:
        raise SystemExit("pip install astroquery astropy to query public MAST products") from exc
    from astropy.time import Time
    start = Time("2026-11-04T00:00:00").mjd
    stop = Time("2026-11-19T00:00:00").mjd
    ids = sorted(set(re.sub(r"^\D+-", "", r["program"]) for r in SOURCES["mast"]["rows"]))
    product_rows = []
    obs_rows = []
    for program in ids:
        observations = Observations.query_criteria(obs_collection="JWST", proposal_id=program)
        chosen = []
        for obs in observations:
            try:
                release = float(obs["t_release"])
                public = str(obs["dataRights"]).upper() == "PUBLIC"
                inst = "NIRSPEC" in str(obs["instrument_name"]).upper()
                if public and inst and start <= release < stop:
                    chosen.append(obs)
                    obs_rows.append({
                        "program": program, "obsid": str(obs["obsid"]),
                        "t_release_mjd": release, "dataRights": "PUBLIC",
                        "instrument_name": str(obs["instrument_name"])
                    })
            except (KeyError, TypeError, ValueError):
                continue
        if chosen:
            from astropy.table import Table
            products = Observations.get_product_list(Table(rows=chosen))
            for p in products:
                name = str(p["productFilename"])
                uri = str(p["dataURI"])
                if (name.lower().endswith(("_x1d.fits", "_s2d.fits", "_cal.fits", "_rate.fits"))
                        and uri.startswith("mast:")):
                    product_rows.append({
                        "program": program, "obsid": str(p["obsID"]),
                        "filename": name, "dataURI": uri,
                        "dataRights": str(p["dataRights"]) if "dataRights" in p.colnames else "OBSERVATION_PUBLIC",
                        "local_sha256": "", "retrieved": "false",
                    })
        print(f"{program}: selected public observations={len(chosen)}")
    args.cache.mkdir(parents=True, exist_ok=True)
    for name, rows, columns in [
        ("mast_observations.csv", obs_rows, ("program", "obsid", "t_release_mjd", "dataRights", "instrument_name")),
        ("mast_products.csv", product_rows, ("program", "obsid", "filename", "dataURI", "dataRights", "local_sha256", "retrieved"))
    ]:
        path = args.cache / name
        with open(path, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
    safe_write_json(args.cache / "mast_query_receipt.json", {
        "queried_utc": utc_now(),
        "release_mjd_window": [start, stop],
        "programs": ids, "public_observation_rows": len(obs_rows),
        "product_rows": len(product_rows),
        "product_manifest_sha256": digests(args.cache / "mast_products.csv")["sha256"],
        "limitation": "Product list only; native FITS pixels not downloaded, extracted or fit.",
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("download", "freeze", "mast-index"):
        p = sub.add_parser(command)
        p.add_argument("--cache", type=Path, required=True)
        if command == "download":
            group = p.add_mutually_exclusive_group(required=True)
            group.add_argument("--all", action="store_true")
            group.add_argument("--catalog-only", action="store_true")
    args = parser.parse_args()
    {"download": cmd_download, "freeze": cmd_freeze, "mast-index": cmd_mast_index}[args.command](args)


if __name__ == "__main__":
    main()
