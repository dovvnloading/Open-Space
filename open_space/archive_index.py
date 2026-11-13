"""Retrieve bounded PDS observational or HEASARC pointing inventories with provenance."""

import argparse
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import re
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

from .archive import write_immutable

PDS = "https://pds.nasa.gov/api/search/1/products"
HEASARC = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
EXOPLANETS = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"


def retrieve(url, cache, accept):
    """Cache exact response bytes by request URL. Use a new directory for an updated snapshot."""
    key = hashlib.sha256(url.encode()).hexdigest()
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / (key + ".response")
    if path.exists():
        raw = path.read_bytes()
    else:
        request = Request(
            url, headers={"Accept": accept, "User-Agent": "Open-Space-Research/0.2"}
        )
        with urlopen(request, timeout=60) as stream:
            raw = stream.read(16_000_001)
        if len(raw) > 16_000_000:
            raise ValueError("Response too large for this query. Reduce its scope.")
    return raw, path


def record_response(raw, path, url):
    # Call only after validating structure and service-level success.
    write_immutable(path, raw)
    return {
        "query_url": url,
        "response_sha256": hashlib.sha256(raw).hexdigest(),
        "bytes": len(raw),
        "cached_at_utc": dt.datetime.fromtimestamp(
            path.stat().st_mtime, dt.timezone.utc
        ).isoformat(timespec="seconds"),
    }


def pds(target, limit, cache):
    if not re.fullmatch(r"urn:nasa:pds:context:target:[a-z0-9_.-]+", target):
        raise ValueError("Expected a PDS target LID")
    query = f'((ref_lid_target eq "{target}") and (product_class eq "Product_Observational"))'
    url = PDS + "?" + urlencode({"q": query, "limit": limit})
    raw, path = retrieve(url, cache, "application/json")
    data = json.loads(raw)
    records = [
        {
            "archive": "PDS",
            "id": row["id"],
            "title": row["title"],
            "start": row.get("start_date_time"),
            "stop": row.get("stop_date_time"),
            "label_url": row["metadata"]["label_url"],
            "data_urls": row.get("properties", {}).get(
                "ops:Data_File_Info.ops:file_ref", []
            ),
            "instruments": [
                t["id"] for t in row.get("observing_system_components", [])
            ],
            "targets": [t["id"] for t in row.get("targets", [])],
        }
        for row in data["data"]
    ]
    return {
        "query": query,
        "total_matches_reported": data["summary"]["hits"],
        "observations": records,
        "source": record_response(raw, path, url),
    }


def parse_votable(raw):
    doc = ET.fromstring(raw)
    infos = [e for e in doc.findall(".//{*}INFO") if e.get("name") == "QUERY_STATUS"]
    if not infos or any(e.get("value") != "OK" for e in infos):
        raise ValueError(
            "TAP query failed or was truncated: "
            + " ".join(e.text or "" for e in infos)
        )
    table = doc.find(".//{*}TABLE")
    if table is None:
        raise ValueError("TAP response missing table")
    if table.find(".//{*}TABLEDATA") is None:
        from astropy.io.votable import parse_single_table
        from io import BytesIO

        parsed = parse_single_table(BytesIO(raw))
        array = parsed.to_table(use_names_over_ids=True).as_array()
        return [
            {
                k.lower(): (
                    row[k].decode() if isinstance(row[k], bytes) else str(row[k])
                )
                for k in array.dtype.names
            }
            for row in array
        ]
    names = [e.get("name").lower() for e in table.findall("{*}FIELD")]
    return [
        dict(zip(names, [td.text or "" for td in tr.findall("{*}TD")], strict=True))
        for tr in table.findall(".//{*}TR")
    ]


def heasarc(ra, dec, radius, limit, cache):
    if not all(math.isfinite(x) for x in (ra, dec, radius)) or not (
        0 <= ra < 360 and -90 <= dec <= 90 and 0 < radius <= 3
    ):
        raise ValueError("Invalid cone coordinates or radius (maximum 3 degrees)")
    query = (
        f"SELECT TOP {limit} obsid,name,ra,dec,exposure_a,exposure_b,public_date FROM numaster "
        "WHERE exposure_a > 0 AND 1=CONTAINS(POINT('ICRS',ra,dec),"
        f"CIRCLE('ICRS',{ra},{dec},{radius})) ORDER BY obsid"
    )
    url = (
        HEASARC
        + "?"
        + urlencode({"request": "doQuery", "lang": "ADQL", "query": query})
    )
    raw, path = retrieve(url, cache, "application/x-votable+xml")
    rows = parse_votable(raw)
    observations = [
        {
            "archive": "HEASARC",
            "mission": "NuSTAR",
            "obsid": row["obsid"],
            "pointing_name": row["name"],
            "ra_deg": float(row["ra"]),
            "dec_deg": float(row["dec"]),
            "module_a_exposure_seconds": float(row["exposure_a"]),
            "module_b_exposure_seconds": float(row["exposure_b"]),
            "public_date_mjd": row["public_date"],
        }
        for row in rows
    ]
    return {
        "query": query,
        "observations": observations,
        "source": record_response(raw, path, url),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, default=Path("runs/inventory-cache"))
    p.add_argument("--limit", type=int, default=5)
    sub = p.add_subparsers(dest="service", required=True)
    planets = sub.add_parser("pds")
    planets.add_argument(
        "target",
        help="PDS logical ID, e.g. urn:nasa:pds:context:target:asteroid.101955_bennu",
    )
    high_energy = sub.add_parser("heasarc")
    high_energy.add_argument("ra", type=float)
    high_energy.add_argument("dec", type=float)
    high_energy.add_argument("--radius", type=float, default=0.25)
    args = p.parse_args()
    if not 1 <= args.limit <= 100:
        p.error("limit must be between 1 and 100")
    try:
        data = (
            pds(args.target, args.limit, args.cache)
            if args.service == "pds"
            else heasarc(args.ra, args.dec, args.radius, args.limit, args.cache)
        )
    except (OSError, KeyError, ValueError, ET.ParseError) as exc:
        p.exit(1, f"inventory: {exc}\n")
    print(json.dumps(data, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
