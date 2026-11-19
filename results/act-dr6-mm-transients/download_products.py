#!/usr/bin/env python3
"""Bounded, checksum-recorded downloads of *exact* frozen MapCat map+ivar pairs.

Requires a previously frozen manifest SHA-256. Does not search images.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
from urllib.request import Request, urlopen

from freeze_mapcat import ROOT, digest


def product_url(relative: str) -> str:
    if not relative or relative.startswith(("http:", "https:")):
        raise ValueError("MapCat must contain bucket-relative paths")
    p = PurePosixPath(relative.lstrip("/"))
    if ".." in p.parts or "?" in relative or "#" in relative:
        raise ValueError("Untrusted MapCat product path")
    return ROOT + "/" + str(p)


def get_fits(url: str, destination: Path, max_bytes: int) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.with_name(destination.name + ".partial")
    if part.exists():
        part.unlink()
    request = Request(url, headers={"User-Agent": "Open-Space-ACT-receipts/1.0"})
    try:
        with urlopen(request, timeout=120) as resp, part.open("wb") as output:
            length = int(resp.headers.get("Content-Length") or "0")
            if length > max_bytes:
                raise ValueError(f"FITS product larger than bound {max_bytes} bytes")
            h = hashlib.sha256()
            size = 0
            while chunk := resp.read(1024 * 1024):
                size += len(chunk)
                if size > max_bytes:
                    raise ValueError("FITS product exceeded download bound")
                h.update(chunk)
                output.write(chunk)
            etag = resp.headers.get("ETag")
            lastmod = resp.headers.get("Last-Modified")
        with part.open("rb") as stream:
            if stream.read(9) != b"SIMPLE  =":
                raise ValueError("Not a primary FITS image; no verified observation")
        part.replace(destination)
        return {"url": url, "sha256": h.hexdigest(), "bytes": size,
                "etag_untrusted": etag, "last_modified_http": lastmod}
    except BaseException:
        part.unlink(missing_ok=True)
        raise


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--manifest-sha256", required=True, help="SHA-256 from frozen inventory receipt")
    ap.add_argument("--map-name", required=True, help="Exact MapCat map_name; no unfrozen discovery")
    ap.add_argument("--out", type=Path, default=Path("act-fits-stage"))
    ap.add_argument("--max-product-mib", type=int, default=512)
    args = ap.parse_args()
    expected = args.manifest_sha256.lower()
    if len(expected) != 64 or digest(args.manifest) != expected:
        ap.error("manifest SHA-256 mismatch; do not silently update the survey selection")
    if not 1 <= args.max_product_mib <= 4096:
        ap.error("max-product-mib must be between 1 and 4096")
    selected = []
    with args.manifest.open(encoding="utf-8") as stream:
        for line in stream:
            row = json.loads(line)
            if row["map_name"] == args.map_name:
                selected.append(row)
    if len(selected) != 1:
        ap.error(f"Expected exactly one frozen map_name; found {len(selected)}")
    row = selected[0]
    products = {}
    for product, column in (("map", "map_path"), ("ivar", "ivar_path")):
        if not row.get(column):
            raise ValueError(f"Catalog row missing mandatory {column}")
        url = product_url(row[column])
        local = args.out / f"{row['map_id']}_{product}.fits"
        products[product] = get_fits(url, local, args.max_product_mib * 1024 * 1024)
        products[product]["local_path"] = str(local)
    receipt = {
        "status": "fits_pair_downloaded_not_searched",
        "frozen_map_name": args.map_name,
        "map_id": str(row["map_id"]),
        "manifest_sha256": expected,
        "utc": datetime.now(timezone.utc).isoformat(),
        "products": products,
        "reference_coadd_staged": False,
        "positive_control_recovery": False,
        "rates_measured": False,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    destination = args.out / f"{row['map_id']}_download_receipt.json"
    destination.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
