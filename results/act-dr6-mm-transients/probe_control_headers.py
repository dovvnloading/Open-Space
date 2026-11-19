#!/usr/bin/env python3
"""Inspect bounded public FITS headers near the three Naess et al. controls.

This identifies which frozen depth-1 products *might cover* published controls.
It does NOT read science pixels or recover a transient; HTTP range bytes are not
full-file SHA receipts. Failures remain in the JSON record.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

from astropy.io import fits
from astropy.wcs import WCS

from freeze_mapcat import ROOT, digest

EXPECTED_SHA = "b6a0c020e3b03bc23aac66536ae1b31cb8889cbc3295c1d2acf16f8de3952503"
# (right ascension, declination) in ICRS degrees, date of peak in the abstract.
CONTROLS = [
    {"name": "Naess2021_2018", "date": "2018-09-11", "ra_deg": 301.9952, "dec_deg": 16.1652},
    {"name": "Naess2021_2019_Nov", "date": "2019-11-08", "ra_deg": 273.8138, "dec_deg": -49.4628},
    {"name": "Naess2021_2019_Dec", "date": "2019-12-15", "ra_deg": 105.1584, "dec_deg": -11.2434},
]
MAX_HEADER_BYTES = 65536


def header_for(url: str) -> tuple[fits.Header, dict]:
    req = Request(url, headers={"Range": f"bytes=0-{MAX_HEADER_BYTES - 1}",
                                "User-Agent": "Open-Space-ACT-positive-control-header-probe/1.0"})
    with urlopen(req, timeout=35) as resp:
        data = resp.read(MAX_HEADER_BYTES)
        meta = {
            "http_status": int(resp.status),
            "content_range": resp.headers.get("Content-Range"),
            "etag_untrusted": resp.headers.get("ETag"),
            "range_bytes_sha256": hashlib.sha256(data).hexdigest(),
            "bytes_inspected": len(data),
        }
    if not data.startswith(b"SIMPLE  ="):
        raise ValueError("Not a primary FITS HDU header")
    if len(data) < 2880:
        raise ValueError("Incomplete FITS primary header")
    # A FITS card is exactly 80 bytes. Stop at END and pad for strict parser.
    cards = [data[i:i+80] for i in range(0, len(data)//80*80, 80)]
    end = next((i for i, x in enumerate(cards) if x[:8].decode("ascii", "ignore").strip() == "END"), None)
    if end is None:
        raise ValueError("FITS header exceeded header-probe size")
    cards = cards[:end+1]
    header = fits.Header.fromstring(b"".join(cards).decode("ascii"), sep="")
    return header, meta


def probe(task: tuple) -> dict:
    control, row = task
    result = {
        "control": control["name"], "catalog_map_name": row["map_name"],
        "map_id": row["map_id"], "catalog_start_time": row["start_time"],
        "array": row["tube_slot"], "frequency": row["frequency"],
        "header_is_not_full_fits_checksum": True,
    }
    try:
        if (not row["map_path"] or "://" in row["map_path"] or ".." in row["map_path"]):
            raise ValueError("Invalid catalog map path")
        url = ROOT + "/" + row["map_path"].lstrip("/")
        h, meta = header_for(url)
        result.update(meta)
        celestial = WCS(h).celestial
        # CRPIX in FITS is one-based; origin 0 returns zero-based pixels.
        x, y = celestial.all_world2pix(control["ra_deg"], control["dec_deg"], 0)
        nx, ny = int(h.get("NAXIS1", 0)), int(h.get("NAXIS2", 0))
        result.update({
            "nx": nx, "ny": ny, "x_pix": round(float(x), 3), "y_pix": round(float(y), 3),
            "within_fits_bounding_box": bool(0 <= x < nx and 0 <= y < ny),
            "wcs_header_validated": True,
        })
    except Exception as e:
        result["error"] = f"{type(e).__name__}: {str(e)[:220]}"
        result["within_fits_bounding_box"] = False
        result["wcs_header_validated"] = False
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--max-maps-per-control", type=int, default=60)
    a = p.parse_args()
    if digest(a.manifest) != EXPECTED_SHA:
        p.error("Eligible manifest does not match pinned SHA-256")
    if not (1 <= a.workers <= 12 and 1 <= a.max_maps_per_control <= 100):
        p.error("Unbounded concurrency/selection")
    rows = [json.loads(s) for s in gzip.open(a.manifest, "rt")]
    targets = []
    eligible = []
    for c in CONTROLS:
        t = datetime.fromisoformat(c["date"])
        possible = [row for row in rows if row["frequency"] == "f150"
                    and abs((datetime.fromisoformat(row["start_time"]) - t).total_seconds()) < 2*86400]
        possible.sort(key=lambda row: (abs((datetime.fromisoformat(row["start_time"]) - t).total_seconds()),
                                       row["map_name"]))
        selection = possible[:a.max_maps_per_control]
        eligible.append({"name": c["name"], "within_two_days_f150_catalog_rows": len(possible),
                         "header_probed_rows": len(selection)})
        targets.extend((c, r) for r in selection)
    observations = []
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures = [pool.submit(probe, target) for target in targets]
        for f in as_completed(futures):
            observations.append(f.result())
    observations.sort(key=lambda r: (r["control"], r["catalog_map_name"]))
    summary = {
        "classification": "published_positive_control_FITS_HEADER_probe_not_science_pixels",
        "source": "Naess et al. 2021, DOI 10.3847/1538-4357/abfe6d",
        "eligible_manifest_sha256": EXPECTED_SHA,
        "controls": eligible,
        "requests": len(observations),
        "valid_fits_wcs_headers": sum(v["wcs_header_validated"] for v in observations),
        "bounding_box_hits": sum(v["within_fits_bounding_box"] for v in observations),
        "errors": sum("error" in v for v in observations),
        "observations": observations,
        "full_fits_bytes_downloaded": False,
        "transient_recovered": False,
        "scientific_result": False,
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "observations"}, indent=2, sort_keys=True))
    print("hits", [dict(control=v["control"], map_name=v["catalog_map_name"])
                   for v in observations if v["within_fits_bounding_box"]][:20])


if __name__ == "__main__":
    main()
