"""Fail-closed provenance and census for the frozen TESS screen.

Pure standard-library functions so cache and ledger guards can be tested offline.
A complete screening census is not a completed scientific investigation.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

PRODUCT = re.compile(r"tess\d+-s(\d{4})-(\d{16})-\d+-s_lc\.fits")
DIGEST = re.compile(r"[a-f0-9]{64}")


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return digest(json.dumps(value, sort_keys=True, separators=(",", ":"),
                             allow_nan=False).encode("utf-8"))


def checked_catalog(path):
    index = json.loads(Path(path).read_bytes())
    if (index.get("status") != "complete"
            or index.get("required_failures") != []
            or index.get("parse_issues") != []):
        raise ValueError("Incomplete catalog exclusion index; refusing to freeze a plan")
    ids = index.get("tic_exclusions")
    if (not isinstance(ids, list) or not ids
            or any(type(tic) is not int or tic <= 0 for tic in ids)
            or len(set(ids)) != len(ids)
            or index.get("unique_tic_exclusion_count") != len(ids)):
        raise ValueError("Malformed or empty exclusion index")
    return index


def product_identity(row):
    tic, sector = row["tic"], row["sector"]
    match = PRODUCT.fullmatch(row["filename"])
    if (type(tic) is not int or type(sector) is not int or tic <= 0
            or not match or int(match[1]) != sector or int(match[2]) != tic):
        raise ValueError("Invalid frozen product identity")
    return tic, sector


def record_name(row):
    tic, sector = product_identity(row)
    return f"tic{tic}_s{sector}.json"


def load_plan(directory):
    directory = Path(directory)
    raw = (directory / "plan.json").read_bytes()
    plan = json.loads(raw)
    targets = plan["targets"]
    if not targets or plan["eligible_targets"] != len(targets):
        raise ValueError("Frozen eligible denominator mismatch")
    tics = set()
    for rank, row in enumerate(targets):
        tic, sector = product_identity(row)
        if (tic in tics or sector != plan["sector"] or row.get("rank") != rank
                or row.get("selection") != "eligible"):
            raise ValueError("Duplicate, unranked or ineligible frozen target")
        tics.add(tic)
    if targets != sorted(targets, key=lambda r: (
            digest(f"open-space-20261125:{r['tic']}".encode()), r["tic"])):
        raise ValueError("Frozen hash ordering mismatch")
    by_tic = {r["tic"]: r for r in targets}
    seen, eligible = set(), set()
    with (directory / "selection.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, strict=True)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Malformed selection header")
        for row in reader:
            if None in row or None in row.values():
                raise ValueError("Malformed selection row")
            tic = int(row["tic"])
            if tic in seen or int(row["sector"]) != plan["sector"]:
                raise ValueError("Duplicate or wrong-sector selection row")
            product_identity({"tic": tic, "sector": int(row["sector"]),
                              "filename": row["filename"]})
            seen.add(tic)
            if row["selection"] == "eligible":
                target = by_tic.get(tic)
                if target is None or any(
                    row.get(k) != ("" if v is None else str(v))
                    for k, v in target.items()
                ):
                    raise ValueError("Eligible row differs from frozen plan")
                eligible.add(tic)
            elif tic in by_tic:
                raise ValueError("Frozen target marked excluded")
    if (eligible != tics or len(seen) != plan["inventory_targets"]
            or plan["inventory"]["n"] != len(seen)):
        raise ValueError("Selection ledger denominator mismatch")
    return plan, digest(raw)


def checked_product(row, meta):
    tic, sector = product_identity(row)
    if (meta.get("tic") != tic or meta.get("sector") != sector
            or meta.get("filename") != row["filename"]
            or meta.get("mast_uri") != "mast:TESS/product/" + row["filename"]
            or type(meta.get("bytes")) is not int or meta["bytes"] <= 0
            or not all(isinstance(meta.get(k), str) and DIGEST.fullmatch(meta[k])
                       for k in ("sha256", "columns_sha256"))):
        raise ValueError("Missing or inconsistent photometry receipt")
    return {k: meta[k] for k in ("filename", "mast_uri", "sha256", "columns_sha256", "bytes")}


def binding(row, plan_sha, configuration_sha, product=None):
    return {"plan_sha256": plan_sha, "configuration_sha256": configuration_sha,
            "target_sha256": canonical(row), "product": product}


def check_record(row, record, plan_sha, configuration_sha, expected_product=None):
    """Raise on stale, foreign, malformed, or partly failed records."""
    tic, sector = product_identity(row)
    if (record.get("tic") != tic or record.get("sector") != sector
            or record.get("rank") != row["rank"]):
        raise ValueError("Screening identity mismatch")
    provenance = record.get("provenance")
    if not isinstance(provenance, dict):
        raise ValueError("Unbound legacy screening record")
    product = provenance.get("product")
    status = record.get("status")
    if product is not None:
        # Validate even during artifact-only aggregation.
        verified = checked_product(row, {**product, "tic": tic, "sector": sector})
        if product != verified:
            raise ValueError("Malformed screening product binding")
    if provenance != binding(row, plan_sha, configuration_sha, product):
        raise ValueError("Plan, target or configuration mismatch")
    if expected_product is not None and product != expected_product:
        raise ValueError("Cached photometry changed")
    peaks = record.get("peaks")
    if not isinstance(peaks, list):
        raise ValueError("Missing peak list")
    if status == "searched":
        if product is None or len(peaks) != 5:
            raise ValueError("Incomplete successful screen")
        for rank, peak in enumerate(peaks, 1):
            if (not isinstance(peak, dict) or peak.get("peak_rank") != rank
                    or type(peak.get("preliminary_pass")) is not bool
                    or any(type(peak.get(k)) not in (int, float)
                           for k in ("period_days", "epoch_btjd", "duration_days", "white_snr"))
                    or not 0 < peak["duration_days"] < peak["period_days"]):
                raise ValueError("Malformed screening peak")
            canonical(peak)  # disallow NaN and infinity
    elif status == "unsearchable":
        if product is None or peaks or not record.get("reason"):
            raise ValueError("Invalid coverage exclusion")
    elif status == "error":
        if peaks or not record.get("error"):
            raise ValueError("Partially failed screen contains usable peaks")
    else:
        raise ValueError("Unknown screening outcome")
    return status


def census(plan, plan_sha, configuration_sha, screen):
    """Count every frozen target once; never promote peaks from a failed record."""
    screen = Path(screen)
    ledger, shortlist = [], []
    expected = {record_name(row) for row in plan["targets"]}
    for row in plan["targets"]:
        name = record_name(row)
        path = screen / name
        item = {"tic": row["tic"], "sector": row["sector"], "rank": row["rank"],
                "filename": row["filename"], "record": name}
        if not path.is_file():
            item["status"] = "missing"
        else:
            try:
                raw = path.read_bytes()
                item["record_sha256"] = digest(raw)
                record = json.loads(raw)
                status = check_record(row, record, plan_sha, configuration_sha)
                item.update(status=status, product=record["provenance"]["product"])
                if status == "searched":
                    passed = [p for p in record["peaks"] if p["preliminary_pass"]]
                    item["preliminary_peaks"] = len(passed)
                    shortlist.extend({"tic": row["tic"], "sector": row["sector"],
                                      "rank": row["rank"], **p} for p in passed)
                elif status == "error":
                    item["reason"] = record["error"]
                else:
                    item["reason"] = record["reason"]
            except (ValueError, TypeError, KeyError, OSError, AttributeError) as exc:
                item.update(status="invalid", reason=str(exc))
        ledger.append(item)
    extras = sorted(p.name for p in screen.glob("*.json") if p.name not in expected)
    counts = Counter(row["status"] for row in ledger)
    shortlist.sort(key=lambda p: (-p["white_snr"], p["rank"], p["peak_rank"]))
    complete = (counts["searched"] + counts["unsearchable"] == len(ledger)
                and not extras)
    summary = {"status": "complete_screening_census" if complete else "incomplete",
               "plan_sha256": plan_sha, "configuration_sha256": configuration_sha,
               "expected_targets": len(ledger),
               "counts": {k: counts[k] for k in
                          ("searched", "unsearchable", "error", "missing", "invalid")},
               "foreign_records": extras, "preliminary_peaks": len(shortlist),
               "preliminary_targets": len({p["tic"] for p in shortlist}),
               "limitations": ["Screening triggers are not candidates or discoveries.",
                               "A complete census is not a completed novelty or false-positive review."]}
    return summary, ledger, shortlist
