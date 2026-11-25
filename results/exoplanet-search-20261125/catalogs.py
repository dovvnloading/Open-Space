"""Freeze candidate catalogs and build a conservative, auditable TIC exclusion index.

Run with Python 3.10+ (standard library for the baseline; Astropy for MRT tables):
    python catalogs.py --output /path/to/catalog-cache

Failure is explicit: a failed download/parse never becomes an empty clean catalog.
The index excludes whole known-candidate stars for a blind first screen; it does
not prohibit a subsequent independently vetted search for additional planets.
Kepler/K2 entries lacking TIC IDs remain available for coordinate/alias checks.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
PINNED_EUMEMIC = "e2eee928fd621f75c1621f6b4206d34bb222ae59"


def stamp():
    return datetime.now(timezone.utc).isoformat()


def tap(query):
    return TAP + "?" + urlencode({"query": query, "format": "csv"})


def fetch(output, name, url, *, query=None, kind="csv", required=True):
    receipt = {"name": name, "url": url, "query": query, "required": required,
               "retrieval_started_utc": stamp(), "status": "pending"}
    for attempt in range(3):
        try:
            req = Request(url, headers={"User-Agent": "Open-Space exoplanet catalog audit/1.0"})
            with urlopen(req, timeout=180) as response:
                payload = response.read()
                receipt.update(http_status=response.status,
                               final_url=response.url,
                               content_type=response.headers.get("Content-Type"),
                               last_modified=response.headers.get("Last-Modified"),
                               etag=response.headers.get("ETag"))
            (output / name).write_bytes(payload)
            receipt.update(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest(),
                           retrieval_finished_utc=stamp())
            if kind == "csv":
                reader = csv.DictReader(io.StringIO(payload.decode("utf-8-sig")))
                rows = list(reader)
                if not reader.fieldnames or len(reader.fieldnames) < 2:
                    raise ValueError("Not a multi-column CSV table")
                if not rows:
                    raise ValueError("Empty catalog is not accepted as a successful snapshot")
                if any(None in row for row in rows):
                    raise ValueError("CSV rows exceed the declared column count")
                receipt.update(row_count=len(rows), columns=reader.fieldnames)
            elif kind == "mrt":
                from astropy.io import ascii
                table = ascii.read(payload.decode("utf-8"), format="cds")
                rows = [{col: "" if str(row[col]) == "--" else str(row[col])
                         for col in table.colnames} for row in table]
                if not rows:
                    raise ValueError("Empty MRT catalog is not accepted as a successful snapshot")
                receipt.update(row_count=len(rows), columns=table.colnames)
            elif kind == "json":
                rows = json.loads(payload)
            else:
                rows = payload.decode("utf-8")
            receipt["status"] = "ok"
            return receipt, rows
        except Exception as exc:
            receipt.update(status="error", error=repr(exc), attempt=attempt + 1,
                           retrieval_finished_utc=stamp())
            if attempt < 2:
                time.sleep(attempt + 1)
    return receipt, None


def parse_tic(value):
    value = str(value or "").strip()
    match = re.fullmatch(r"(?:TIC\s*)?(\d+)(?:\.0)?", value, re.I)
    return int(match.group(1)) if match and int(match.group(1)) > 0 else None


def supplemental(output):
    """Write a separate literature index without modifying baseline snapshots."""
    jobs = [
        ("raven_nsfp09.csv", "https://zenodo.org/records/19661443/files/nsfp09_table.csv?download=1", "csv"),
        ("raven_vet.csv", "https://zenodo.org/records/19661443/files/vet_table.csv?download=1", "csv"),
        ("raven_val.csv", "https://zenodo.org/records/19661443/files/val_table.csv?download=1", "csv"),
        ("raven_val_rp8.csv", "https://zenodo.org/records/19661443/files/val_rp8_table.csv?download=1", "csv"),
        ("ultrasat_candidates.txt", "https://zenodo.org/records/23504163/files/candidate_table.txt?download=1", "mrt"),
        ("ultrasat_tce.txt", "https://zenodo.org/records/23504163/files/tce_table.txt?download=1", "mrt"),
        ("t16_table1.txt", "https://content.cld.iop.org/journals/0067-0049/284/1/19/revision1/apjsae5b6ct1_mrt.txt", "mrt"),
        ("t16_table2.txt", "https://content.cld.iop.org/journals/0067-0049/284/1/19/revision1/apjsae5b6ct2_mrt.txt", "mrt"),
    ]
    receipts, all_rows, source_counts, tic_sources, parse_issues = [], [], {}, {}, []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(fetch, output, name, url, kind=kind) for name, url, kind in jobs]
        for future in as_completed(futures):
            receipt, rows = future.result()
            receipts.append(receipt)
            name = receipt["name"]
            if rows is None:
                print(name, "error", receipt.get("error"), flush=True)
                continue
            field = next((c for c in ["ticid", "TICID", "TIC", "TIC ID"] if c in receipt["columns"]), None)
            if field is None:
                parse_issues.append({"source": name, "error": "No recognized TIC column"})
                continue
            source_tics = set()
            for number, row in enumerate(rows, 1):
                tic = parse_tic(row.get(field))
                if tic is not None:
                    source_tics.add(tic)
                    tic_sources.setdefault(str(tic), []).append(name)
                elif row.get(field, "").strip() not in ("", "--", "nan", "NaN", "0"):
                    parse_issues.append({"source": name, "row": number, "error": "Malformed TIC", "value": row.get(field)})
                all_rows.append({"source": name, "source_row": number, "tic": tic, "data": row})
            source_counts[name] = {"row_count": len(rows), "tic_column": field,
                                   "unique_tic_count": len(source_tics),
                                   "rows_without_tic": sum(parse_tic(row.get(field)) is None for row in rows)}
            print(name, len(rows), len(source_tics), flush=True)
    failures = [r["name"] for r in receipts if r["status"] != "ok"]
    index = {
        "created_utc": stamp(),
        "status": "complete" if not failures and not parse_issues else "incomplete",
        "purpose": "Supplemental published-signal TIC exclusion index; union with the separately frozen baseline index.",
        "required_failures": failures, "parse_issues": parse_issues,
        "source_counts": source_counts,
        "tic_exclusions": sorted(map(int, tic_sources)),
        "unique_tic_exclusion_count": len(tic_sources),
        "tic_sources": {key: sorted(set(value)) for key, value in sorted(tic_sources.items(), key=lambda r: int(r[0]))},
        "limitations": [
            "Whole-star exclusion is conservative; not every RAVEN NSFP>0.9 or ULTRASAT TCE is a vetted planet candidate.",
            "Retains all dispositions to avoid rebranding previously reported false alarms or eclipsing binaries.",
            "TIC-less entries require coordinate/Gaia-ID crossmatching; no-hit in this index does not establish novelty.",
            "This supplements rather than overwrites baseline TOI/CTOI/confirmed/KOI/K2/community snapshots.",
        ],
    }
    (output / "catalog_receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")
    (output / "catalog_rows.json").write_text(json.dumps(all_rows, separators=(",", ":")) + "\n")
    (output / "exclusion_index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(json.dumps({k: index[k] for k in ["status", "required_failures", "parse_issues", "source_counts", "unique_tic_exclusion_count"]}, indent=2))
    return 0 if index["status"] == "complete" else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--supplemental-only", action="store_true",
                    help="Acquire T16/RAVEN/ULTRASAT only; use a separate output directory")
    args = ap.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    if args.supplemental_only:
        return supplemental(output)
    queries = {
        "nasa_pscomppars.csv": "select pl_name,hostname,tic_id,gaia_dr2_id,gaia_dr3_id,ra,dec,pl_orbper,pl_tranmid,pl_trandur,discoverymethod,disc_year from pscomppars",
        "nasa_toi.csv": "select toi,tid,ctoi_alias,tfopwg_disp,ra,dec,pl_orbper,pl_tranmid,pl_trandurh from toi",
        "nasa_cumulative.csv": "select kepid,kepoi_name,kepler_name,koi_disposition,koi_pdisposition,koi_period,koi_time0bk,koi_duration,ra,dec from cumulative",
        "nasa_k2pandc.csv": "select * from k2pandc",
    }
    jobs = [(name, tap(query), {"query": query}) for name, query in queries.items()]
    jobs += [
        ("exofop_tois.csv", "https://exofop.ipac.caltech.edu/tess/json/tois.csv", {}),
        ("exofop_ctois.csv", "https://exofop.ipac.caltech.edu/tess/json/ctois.csv", {}),
        ("eumemic_pinned_ephemerides.csv", f"https://raw.githubusercontent.com/eumemic/exoplanets/{PINNED_EUMEMIC}/candidates/ephemerides.csv", {}),
        ("eumemic_repository.json", "https://api.github.com/repos/eumemic/exoplanets", {"kind": "json"}),
        ("eumemic_head.json", "https://api.github.com/repos/eumemic/exoplanets/commits/HEAD", {"kind": "json"}),
        ("rpinho_tree.json", "https://api.github.com/repos/rpinho/planet-hunt/git/trees/HEAD?recursive=1", {"kind": "json", "required": False}),
    ]
    snapshots = {}
    receipts = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(fetch, output, name, url, **kwargs): name for name, url, kwargs in jobs}
        for future in as_completed(futures):
            receipt, rows = future.result()
            receipts.append(receipt)
            if rows is not None:
                snapshots[receipt["name"]] = rows
            print(receipt["name"], receipt["status"], receipt.get("row_count", receipt.get("bytes")), flush=True)
    head = snapshots.get("eumemic_head.json", {}).get("sha")
    if head:
        receipt, rows = fetch(output, "eumemic_current_ephemerides.csv", f"https://raw.githubusercontent.com/eumemic/exoplanets/{head}/candidates/ephemerides.csv")
        receipt["commit_sha"] = head
        receipts.append(receipt)
        if rows is not None:
            snapshots[receipt["name"]] = rows
        print(receipt["name"], receipt["status"], receipt.get("row_count"), flush=True)
    rpinho_tree = snapshots.get("rpinho_tree.json", {})
    rpinho_head = rpinho_tree.get("sha")
    rpinho_paths = [r["path"] for r in rpinho_tree.get("tree", []) if r.get("type") == "blob" and (r["path"].lower().endswith((".csv", ".json", ".md"))) and int(r.get("size", 0)) < 3000000]
    for path in rpinho_paths:
        name = "rpinho_" + path.replace("/", "__")
        kind = "csv" if path.endswith(".csv") else "json" if path.endswith(".json") else "text"
        receipt, rows = fetch(output, name, f"https://raw.githubusercontent.com/rpinho/planet-hunt/{rpinho_head}/{path}", kind=kind, required=False)
        receipt["commit_sha"] = rpinho_head
        receipts.append(receipt)
        if rows is not None:
            snapshots[name] = rows
    tic_sources = {}
    source_counts = {}
    row_index = []
    parse_issues = []
    for receipt in receipts:
        name = receipt["name"]
        rows = snapshots.get(name)
        if not name.endswith(".csv") or not isinstance(rows, list):
            continue
        fields = receipt["columns"]
        tic_column = next((col for col in ["TIC ID", "TIC", "tic_id", "tid", "tic"] if col in fields), None)
        source_tics = set()
        malformed = []
        for number, row in enumerate(rows, 2):
            value = row.get(tic_column) if tic_column else None
            tic = parse_tic(value)
            if tic is not None:
                source_tics.add(tic)
                tic_sources.setdefault(str(tic), []).append(name)
            elif value and value.strip() not in ("", "null", "NULL", "nan", "NaN", "--"):
                malformed.append({"row": number, "value": value})
            row_index.append({"source": name, "source_row": number, "tic": tic, "data": row})
        source_counts[name] = {"row_count": len(rows), "tic_column": tic_column,
                               "unique_tic_count": len(source_tics),
                               "rows_without_tic": sum(parse_tic(r.get(tic_column)) is None for r in rows),
                               "malformed_tic_count": len(malformed)}
        if malformed:
            parse_issues.append({"source": name, "issue": "nonempty_malformed_tic", "examples": malformed[:10], "count": len(malformed)})
    # Readme/discussion candidates are not a homogeneous catalog; extract only
    # explicitly prefixed TIC identifiers and preserve their source context.
    text_tic_sources = {}
    for name, contents in snapshots.items():
        if name.startswith("rpinho_") and isinstance(contents, str):
            for match in re.finditer(r"\bTIC[\s_-]*(\d{5,12})\b", contents, re.I):
                tic = int(match.group(1))
                tic_sources.setdefault(str(tic), []).append(name)
                text_tic_sources.setdefault(str(tic), []).append({"source": name, "context": contents[max(0, match.start()-80):match.end()+120]})
    required_failures = [r["name"] for r in receipts if r["required"] and r["status"] != "ok"]
    summary = {
        "created_utc": stamp(),
        "purpose": "Conservative whole-TIC exclusion for a blind new-candidate screen; not proof of novelty for nonexcluded stars.",
        "status": "complete" if not required_failures and not parse_issues else "incomplete",
        "required_failures": required_failures,
        "parse_issues": parse_issues,
        "source_counts": source_counts,
        "tic_exclusions": sorted(map(int, tic_sources)),
        "unique_tic_exclusion_count": len(tic_sources),
        "tic_sources": {key: sorted(set(value)) for key, value in sorted(tic_sources.items(), key=lambda r: int(r[0]))},
        "rpinho_text_tic_evidence": text_tic_sources,
        "limitations": [
            "Kepler/K2 entries may lack TIC IDs; inspect retained coordinates and aliases for any promoted candidate.",
            "Not an exhaustive literature index, unpublished work is unknowable, and indexed literature can lag.",
            "A matching star excludes this blind screen, but an additional independent period would need event-level comparison.",
            "All dispositions are retained, including known false positives and false alarms.",
        ],
    }
    (output / "receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")
    (output / "catalog_rows.json").write_text(json.dumps(row_index, separators=(",", ":")) + "\n")
    (output / "exclusion_index.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: summary[key] for key in ["status", "required_failures", "parse_issues", "unique_tic_exclusion_count", "source_counts"]}, indent=2))
    return 0 if summary["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
