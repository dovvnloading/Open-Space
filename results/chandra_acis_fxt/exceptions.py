"""Resolve only missing CXC release dates against the public HEASARC mirror.

Source definitions
https://cda.harvard.edu/cxctap/tables
https://heasarc.gsfc.nasa.gov/W3Browse/all/chanmaster.html
https://heasarc.gsfc.nasa.gov/docs/archive/vo/instructions.html
"""
import argparse
from collections import defaultdict
from datetime import date, timedelta
from decimal import Decimal, InvalidOperation, ROUND_FLOOR
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

from astropy.io.votable import parse_single_table

from archive import download, dump, FIELDS, PROPOSAL_FIELDS, TAP

HEASARC = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
COLUMNS = ("obsid", "public_date", "status", "detector", "grating", "data_mode", "exposure")


def read_json(path):
    return json.loads(path.read_bytes())


def sha(data):
    return hashlib.sha256(data).hexdigest()


def mirror_date(value):
    """Public_Date is normally MJD. Recognize ISO displays and reject sentinels."""
    text = str(value).strip()
    try:
        if re.match(r"^\d{4}-\d{2}-\d{2}(?:$|T| )", text):
            released = date.fromisoformat(text[:10])
        else:
            mjd = Decimal(text)
            if not mjd.is_finite() or not 0 < mjd < 200000:
                return None
            released = date(1858, 11, 17) + timedelta(days=int(mjd.to_integral_value(rounding=ROUND_FLOOR)))
    except (ValueError, InvalidOperation, OverflowError):
        return None
    return released if released >= date(1999, 7, 23) else None


def cxc_metadata(folder, sources, ids):
    """Read the hashed raw census slabs already committed by the survey."""
    result = defaultdict(list)
    for s in sources:
        m = re.fullmatch(r"metadata_(\d+)_(\d+)\.xml\.gz", s["raw_gzip"])
        if not m or s.get("overflow"):
            continue
        if not any(int(m[1]) <= n <= int(m[2]) for n in ids):
            continue
        raw = gzip.decompress((folder / s["raw_gzip"]).read_bytes())
        if sha(raw) != s["sha256"] or len(raw) != s["bytes"]:
            raise ValueError(f'frozen CXC source differs from receipt: {s["raw_gzip"]}')
        table = ET.fromstring(raw).find(".//{*}TABLE")
        fields = [f.attrib["name"].lower().split(".")[-1] for f in table.findall("{*}FIELD")]
        for tr in table.findall(".//{*}TR"):
            r = dict(zip(fields, ((td.text or "").strip() for td in tr), strict=True))
            oid = int(r["obsid"])
            if oid in ids:
                result[oid].append(r)
    if ids - result.keys():
        raise ValueError(f"missing frozen CXC metadata rows for {sorted(ids - result.keys())}")
    return result


def fetch_mirror(ids, out):
    # String constants work with string and numeric ObsID column types on TAP
    # services. Include the zero-padded rendering for archival identifiers.
    values = sorted({f"'{n}'" for n in ids} | {f"'{n:05d}'" for n in ids})
    query = ("SELECT " + ", ".join("c." + col for col in COLUMNS) + " FROM chanmaster c "
             + "WHERE c.obsid IN (" + ",".join(values) + ")")
    body = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "votable", "MAXREC": 10000,
                      "QUERY": query}).encode("ascii")
    raw_path = out / "heasarc_chanmaster.xml"
    receipt_path = out / "heasarc_request.json"
    if raw_path.is_file() and receipt_path.is_file():
        receipt = read_json(receipt_path)
        raw = raw_path.read_bytes()
        if (receipt["post_body"] != body.decode() or sha(raw) != receipt["sha256"]
                or len(raw) != receipt["bytes"]):
            raise ValueError("cached supplement response differs from its request or checksum")
    else:
        receipt = download(HEASARC, raw_path, data=body)
        receipt.update(query=query, post_body=body.decode(), response_file=raw_path.name)
        dump(receipt_path, receipt)
        raw = raw_path.read_bytes()
    xml = ET.fromstring(raw)
    for info in xml.findall(".//{*}INFO"):
        if info.get("name", "").upper() == "QUERY_STATUS" and info.get("value", "").upper() != "OK":
            raise ValueError(f"HEASARC TAP {info.get('value')}: {info.text}")
    # HEASARC can serialize VOTables as BINARY/BINARY2 or TABLEDATA. Retain
    # each field's published name and mask rather than decoding a display table.
    table = parse_single_table(str(raw_path)).to_table(use_names_over_ids=True)
    names = {f.lower().split(".")[-1]: f for f in table.colnames}
    if set(COLUMNS) != set(names):
        raise ValueError(f"unexpected HEASARC fields: {list(names)}")
    by_id = defaultdict(list)
    for r in table:
        mirror_row = {}
        for key, name in names.items():
            val = r[name]
            if bool(getattr(val, "mask", False)):
                val = ""
            if isinstance(val, bytes):
                val = val.decode("ascii")
            mirror_row[key] = str(val).strip()
        oid = int(mirror_row["obsid"])
        if oid not in ids:
            raise ValueError(f"HEASARC returned an unrequested ObsID: {oid}")
        by_id[oid].append(mirror_row)
    return by_id, receipt


def classify(rows, first, last, as_of):
    if not rows:
        return "unknown", None, "missing_chanmaster_row"
    for r in rows:
        if r["status"].upper() != "ARCHIVED":
            return "unknown", None, "unverified_public_status"
        if (not r["detector"].upper().startswith("ACIS") or r["grating"].upper() != "NONE"
                or not r["data_mode"].upper().startswith("T")):
            return "unknown", None, "unverified_or_conflicting_imaging_mode"
    dates = {mirror_date(r["public_date"]) for r in rows}
    if None in dates or len(dates) != 1:
        return "unknown", None, "missing_or_conflicting_public_date"
    day = next(iter(dates))
    if day > as_of:
        return "unknown", None, "public_date_after_mirror_retrieval"
    if first <= day <= last:
        return "eligible", day.isoformat(), "heasarc_archived_release_in_window"
    return "ineligible", day.isoformat(), "heasarc_archived_release_outside_window"


def supplement(inventory_dir, output_dir):
    folder, out = Path(inventory_dir), Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    filenames = ("inventory.json", "inventory_exceptions.json", "inventory_sources.json")
    inputs = {n: sha((folder / n).read_bytes()) for n in filenames}
    existing = {r["obsid"] for r in read_json(folder / filenames[0])}
    exceptions = read_json(folder / filenames[1])
    sources = read_json(folder / filenames[2])
    first = date.fromisoformat(sources["start_inclusive"])
    last = date.fromisoformat(sources["end_inclusive"])
    unresolved = {int(r["obsid"]) for r in exceptions if r["reason"] == "missing_public_date"}
    result = {"eligible": [], "ineligible": [], "unknown": [], "already_covered": sorted(existing & unresolved),
              "frozen_input_sha256": inputs, "start_inclusive": first.isoformat(), "end_inclusive": last.isoformat(),
              "source": "HEASARC CHANMASTER is a twice-weekly CXC mirror. Only a missing CXC release date is resolved."}
    needed = unresolved - existing
    if needed:
        primary = cxc_metadata(folder, sources["sources"], needed)
        mirror, receipt = fetch_mirror(needed, out)
        result["mirror_receipt"] = receipt
        as_of = date.fromisoformat(receipt["retrieved_utc"][:10])
        for oid in sorted(needed):
            rows = [r for r in primary[oid] if r["status"].upper() == "ARCHIVED"
                    and r["grating"].upper() not in {"HETG", "LETG"}
                    and r["exposure_mode"].upper() not in {"CC", "CONTINUOUS", "CONTINUOUS CLOCKING"}]
            catalog = mirror[oid]
            base = {"obsid": oid, "mirror_rows": catalog, "release_date_source": "HEASARC chanmaster.public_date"}
            if not rows or any(r["public_avail_date"] for r in rows):
                result["unknown"].append({**base, "reason": "frozen_cxc_eligibility_requires_review"})
                continue
            status, released, reason = classify(catalog, first, last, as_of)
            base.update(release_date=released, reason=reason)
            if status != "eligible":
                result[status].append(base)
                continue
            row = min(rows, key=lambda r: tuple(r[k] for k in FIELDS + PROPOSAL_FIELDS))
            proposal = row["proposal_number"]
            params = {"REQUEST": "doQuery", "FORMAT": "text", "LANG": "ADQL",
                      "QUERY": "SELECT * FROM cxc.proposal p WHERE p.number = '" + proposal.replace("'", "''") + "'"}
            result["eligible"].append({
                **base, **{k: row[k] for k in ("instrument", "grating", "exposure_mode", "start_date",
                                               "target_name", "proposal_number", *PROPOSAL_FIELDS)},
                "status": "ARCHIVED", "known_release_dates": [released], "metadata_rows": len(rows),
                "event_count": int(row["event_count"]) if row["event_count"] else None,
                "exposure_s": 1000 * float(row["exposure_time"]) if row["exposure_time"] else None,
                "catalog_ra_deg": float(row["ra"]) if row["ra"] else None,
                "catalog_dec_deg": float(row["dec"]) if row["dec"] else None,
                "observation_url": f"https://cda.harvard.edu/chaser/startViewer.do?menuItem=details&obsid={oid}",
                "proposal_url": TAP + "?" + urlencode(params) if proposal else None,
            })
    dump(out / "supplement.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    result = supplement(args.inventory_dir, args.output_dir)
    print("Release-date supplement", {k: len(result[k]) for k in ("eligible", "ineligible", "unknown", "already_covered")})
