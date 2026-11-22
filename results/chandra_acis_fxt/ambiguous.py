"""Search separately the missing-release-date ObsIDs from the frozen survey."""
import argparse
from datetime import date
import hashlib
from html import unescape
import os
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory
import time
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

import campaign
import report
from archive import download, dump
from exceptions import cxc_metadata
from search import DURATIONS

CSC_TAP = "https://cda.cfa.harvard.edu/csc22tap/sync"
CSC_SELECTION = "https://cxc.cfa.harvard.edu/csc2.2/"
CSC_HISTORY = "https://cxc.cfa.harvard.edu/csc2.2/versions.html"
CSC_RELEASE = "2026-10-05"
DARK = re.compile(r"\b(stowed|dark (?:calibration|current|frames?|exposures?|measurement)|bias frames?)\b", re.I)
FROZEN_POPULATION = "172b673a39b7c2bca8b82294563a62a7cb5807fdcc393a8da27d31e56f3c518e"


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def archived_http(folder, label, url, body=None):
    target, receipt = folder / label, folder / (label + ".receipt.json")
    request = body.decode() if body else None
    if target.is_file() and receipt.is_file():
        record = campaign.read(receipt)
        if (record["sha256"] != sha(target) or record.get("post_body") != request or record["url"] != url):
            raise ValueError(f"Altered archived response {label}")
    else:
        record = download(url, target, data=body)
        record.update(file=label, post_body=request)
        dump(receipt, record)
    return target.read_bytes(), record


def frozen_csc(ids, folder):
    """Only a positive per-observation detection in the frozen release counts."""
    folder.mkdir(parents=True, exist_ok=True)
    coverage, sources = {}, {}
    try:
        for label, url in (("selection.html", CSC_SELECTION), ("release_history.html", CSC_HISTORY)):
            raw, sources[label] = archived_http(folder, label, url)
            plain = " ".join(unescape(re.sub(r"<[^>]+>", " ", raw.decode("utf-8", "replace"))).split())
            if label == "selection.html":
                if not (re.search(r"ACIS", plain, re.I) and
                        re.search(r"observations released publicly through the end of 20\d\d", plain, re.I)):
                    raise ValueError("Missing mission statement that CSC 2.2 selects publicly released observations")
            elif not re.search(r"CSC 2\.2\s+2026-10-05", plain, re.I):
                raise ValueError("Cannot verify the frozen CSC 2.2 publication date")
        members = ",".join(map(str, sorted(ids)))
        # The maximum broad-band significance is for an individual observation,
        # not a stack. A conservative >=3 detection supplies positive evidence.
        query = ("SELECT o.obsid, MAX(o.flux_significance_b) AS significance "
                 "FROM csc22.observation_source o WHERE o.instrument='ACIS' "
                 f"AND o.grating='NONE' AND o.obsid IN ({members}) GROUP BY o.obsid")
        body = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "votable",
                          "QUERY": query, "MAXREC": 1000}).encode()
        raw, sources["matches.xml"] = archived_http(folder, "matches.xml", CSC_TAP, body)
        doc = ET.fromstring(raw)
        for info in doc.findall(".//{*}INFO"):
            if info.get("name", "").upper() == "QUERY_STATUS" and info.get("value", "").upper() != "OK":
                raise ValueError(f"CSC TAP {info.get('value')}: {info.text}")
        table = doc.find(".//{*}TABLE")
        if table is None or table.find(".//{*}TABLEDATA") is None:
            raise ValueError("CSC TAP response lacks TABLEDATA")
        fields = [f.attrib["name"].lower().split(".")[-1] for f in table.findall("{*}FIELD")]
        if set(fields) != {"obsid", "significance"}:
            raise ValueError(f"Unexpected CSC columns: {fields}")
        for tr in table.findall(".//{*}TR"):
            row = dict(zip(fields, ((v.text or "").strip() for v in tr), strict=True))
            oid = int(row["obsid"])
            if oid not in ids or oid in coverage:
                raise ValueError("Invalid ObsID in grouped CSC response")
            val = float(row["significance"]) if row["significance"] else None
            coverage[oid] = val if val is not None and -1e100 < val < 1e100 else None
        status, error = "complete", None
    except Exception as exc:
        status, error, coverage = "unknown", f"{type(exc).__name__}: {exc}", {}
    record = {"status": status, "error": error, "sources": sources,
              "release": CSC_RELEASE, "maximum_broad_flux_significance": coverage}
    dump(folder / "lookup.json", record)
    return coverage


def prepare(evidence, output):
    inv = evidence / "inventory"
    supplement_path = evidence / "release_supplement/supplement.json"
    supplement = campaign.read(supplement_path)
    inputs = supplement["frozen_input_sha256"]
    if sha(inv / "inventory.json") != FROZEN_POPULATION:
        raise ValueError("This pass belongs to the original frozen release inventory")
    for name in ("inventory.json", "inventory_exceptions.json", "inventory_sources.json"):
        if sha(inv / name) != inputs[name]:
            raise ValueError(f"Original census input no longer matches {name}")
    freeze = campaign.read(inv / "campaign.json")
    if (freeze["implementation"] != campaign.identity() or freeze["n_controls"] != campaign.N_CONTROL
            or freeze["n_injections_per_product"] != campaign.N_INJECT):
        raise ValueError("The supplementary science implementation differs from the main run")
    ids = {int(item["obsid"]) for item in supplement["unknown"]}
    if ids & {r["obsid"] for r in campaign.read(inv / "inventory.json")}:
        raise ValueError("Supplement would overlap the original dated census")
    lower = date.fromisoformat(supplement["start_inclusive"])
    upper = date.fromisoformat(supplement["end_inclusive"])
    sources = campaign.read(inv / "inventory_sources.json")["sources"]
    primary = cxc_metadata(inv, sources, ids)
    csc = frozen_csc(ids, output / "csc22") if ids else {}
    records = []
    for oid in sorted(ids):
        rows = primary[oid]
        # The earliest start of any archived imaging row can bound acquisition
        # of the ObsID; incomplete starts cannot establish a lower bound.
        archived = [r for r in rows if r["status"].upper() == "ARCHIVED"
                    and r["grating"].upper() not in {"LETG", "HETG"}
                    and r["exposure_mode"].upper() not in {"CC", "CONTINUOUS", "CONTINUOUS CLOCKING"}]
        lower_date = (min(date.fromisoformat(r["start_date"][:10]) for r in archived)
                      if archived and all(r["start_date"] for r in archived) else None)
        significance = csc.get(oid)
        bound = date.fromisoformat(CSC_RELEASE) if significance is not None and significance >= 3 else None
        is_member = bool(lower_date and bound and lower <= lower_date <= bound <= upper)
        targets = sorted({r["target_name"] for r in archived})
        proposals = sorted({r["proposal_type"] for r in archived})
        title = " ".join(targets + [r["proposal_title"] for r in archived]).lower()
        stowed = bool(archived) and all("stowed" in r["target_name"].lower() for r in archived)
        flagged = stowed or bool(DARK.search(title)) or any(
            t.upper() == "CAL" for t in proposals) and bool(re.search(r"\bdark\b", title))
        records.append({"obsid": oid, "target_name": " | ".join(targets), "proposal_type": " | ".join(proposals),
                        "proposal_titles": sorted({r["proposal_title"] for r in archived}), "release_date": None,
                        "earliest_observation_utc_date": str(lower_date) if lower_date else None,
                        "release_upper_bound_utc_date": str(bound) if bound else None,
                        "csc22_max_flux_significance_b": significance, "proven_in_window": is_member,
                        "sky_status": "archived_stowed_target" if stowed else "dark_calibration_review" if flagged
                                      else "nominal_imaging", "raw_cxc_rows": rows})
    dump(output / "release_bounds.json", {"obsids": records, "supplement_input_sha256": sha(supplement_path),
         "implementation": campaign.identity(), "supplemental_program_sha256": sha(Path(__file__)),
         "frozen_catalog_uses_public_observations": CSC_SELECTION, "catalog_release_date_reference": CSC_HISTORY,
         "run_id": os.environ["GITHUB_RUN_ID"], "git_commit": os.environ["GITHUB_SHA"],
         "method": "Earlier than the pinned CSC 2.2 publication when that catalog detects this specific "
                   "ObsID and its selection documentation explicitly requires public observations. "
                   "The observing date supplies only the lower bound. The current revision's historical "
                   "availability and exact release date are unknown. Missing detections supply no bound."})
    return records


def scan_one(oid, output):
    # Give the existing worker the separate prior ledger on retries. This leaves
    # all thresholds, deterministic hash seeds and control/injection draws intact.
    campaign.EVIDENCE = output
    catalog_cache = output / "catalogs" / (os.environ["GITHUB_RUN_ID"] + "-" + os.environ["GITHUB_RUN_ATTEMPT"])
    catalog_cache.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory() as work:
        campaign.process(oid, output, Path(work), catalog_cache)


def reduce_scope(output, name, entries):
    directory = output / name
    directory.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory() as tmp:
        view = Path(tmp)
        (view / "inventory").mkdir()
        (view / "observations").symlink_to(output / "observations", target_is_directory=True)
        dump(view / "inventory/inventory.json", entries)
        dump(view / "inventory/inventory_exceptions.json", [])
        report.EVIDENCE = view
        _, _, observations, trials, tally, exposure, receipts, injections, complete, scanned = report.collect()
    report.rows_to_csv(directory / "coverage.csv", ["obsid", "status", "scanned_all_products", "n_products", "errors"], observations)
    report.rows_to_csv(directory / "event_inputs.csv", ["obsid", "filename", "url", "bytes", "sha256", "fits_sha256",
                       "revision", "included_in_inference"], receipts)
    controls = [report.distinct(report.kept(trials[i])) for i in range(campaign.N_CONTROL)]
    dump(directory / "control_counts.json", [{"replica": i, "retained_or_unresolved": len(c)} for i, c in enumerate(controls)])
    grid = []
    for (duration, photons, psf), b in sorted(tally.items()):
        a, k = b["attempted"], b["recovered"]
        grid.append({"duration_s": duration, "injected_counts": photons, "sigma_sky_px": psf, "planned": b["planned"],
                     "attempted": a, "recovered": k, "fraction": k/a if a else None, "binomial_95": report.bounds(k, a)})
    dump(directory / "injections.json", grid)
    complete_vetting = complete == scanned
    maxima = [max((e["seed_score"] for e in group), default=None) for group in controls]
    if not complete_vetting:
        maxima = [max((e["seed_score"] for e in trials[i] if e["instrument_pass"]), default=None)
                  for i in range(campaign.N_CONTROL)]
    real = report.distinct(report.kept(trials[-1]))
    for lead in real:
        count = sum(val is not None and val >= lead["seed_score"] for val in maxima)
        lead["supplement_fap"] = (1 + count) / (1 + campaign.N_CONTROL)
        lead["reference"] = "retained_or_unresolved" if complete_vetting else "instrument_pass"
    dump(directory / "leads.json", real)
    return {"obsids": [r["obsid"] for r in entries], "fully_scanned": scanned, "followup_complete": complete,
            "eligible_start_support_s": {str(d): exposure[str(d)] for d in DURATIONS},
            "real_raw_triggers": len(trials[-1]), "surviving_or_unresolved": len(real),
            "mean_control_retained_or_unresolved": sum(map(len, controls))/len(controls) if scanned else None,
            "attempted_injections": sum(j["attempted"] for j in injections),
            "recovered_injections": sum(bool(j["recovered"]) for j in injections)}


def summarize(output, entries):
    groups = {"date_bracketed_sky": [], "unresolved_date_sky": [], "non_sky_or_calibration_review": []}
    statuses = []
    for entry in entries:
        oid = entry["obsid"]
        path = output / "observations" / f"{oid}.json.gz"
        ledger = campaign.read_obs(path) if path.is_file() else {"status": "unattempted", "products": [], "errors": []}
        scans = [r["search"] for r in ledger["products"] if r.get("complete") and "search" in r]
        if entry["sky_status"] != "nominal_imaging":
            group = "non_sky_or_calibration_review"
        else:
            group = "date_bracketed_sky" if entry["proven_in_window"] else "unresolved_date_sky"
        groups[group].append(entry)
        full = ledger["status"] in {"complete", "followup_incomplete"} and bool(scans)
        # Do not turn a server/format failure into an apparent observation of zero seconds.
        statuses.append({"obsid": oid, "date_sky_group": group, "worker_status": ledger["status"],
                         "confirmed_zero_60s_support": full and not any(s["exposure_by_duration_s"]["60"] > 0 for s in scans),
                         "usable_products": len(scans), "errors": ledger["errors"]})
    results = {name: reduce_scope(output, name, rows) for name, rows in groups.items()}
    dump(output / "summary.json", {"groups": results, "status": statuses,
         "scope": "Supplemental cohorts only. No time or candidates are included in the dated survey. "
                  "For non-sky or dark-calibration flagged data, numerical start support is a detector control "
                  "only and carries zero claimed celestial survey exposure."})
    lines = ["# Chandra supplemental search for observations with no catalog release date", "",
             f"**{len(entries)} ambiguous ObsIDs were assigned; {sum(r['fully_scanned'] for r in results.values())} "
             "completed the blind scan and controls.**", "",
             "| Cohort | Assigned | Scanned | 60 s start support (s) | Raw triggers | Retained or unresolved |",
             "|---|---:|---:|---:|---:|---:|"]
    for name, r in results.items():
        time_result = r["eligible_start_support_s"]["60"] if name != "non_sky_or_calibration_review" else 0.
        lines.append(f"| [{name}]({name}/leads.json) | {len(r['obsids'])} | {r['fully_scanned']} | "
                     f"{time_result:.1f} | {r['real_raw_triggers']} | {r['surviving_or_unresolved']} |")
    lines.extend(["", "These cohorts are outside the primary dated search. Each cohort records exact event URLs "
                  "and SHA-256 hashes, all 24 scrambled catalogs per product, and a 12-injection-per-product recovery grid. "
                  "The dark-calibration cohort contributes zero claimed sky exposure. Acquisition failures and "
                  "completed scans with zero time support are recorded individually.",
                  "", "[Release bounds](release_bounds.json) · [Coverage and errors](summary.json)"])
    (output / "README.md").write_text("\n".join(lines) + "\n")
    return statuses


def main(args):
    if (os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_REF") != "refs/heads/research/chandra-acis-fast-transients"
            or os.getenv("GITHUB_EVENT_NAME") not in {"push", "workflow_dispatch"}):
        raise RuntimeError("Supplemental acquisition requires the research branch workflow")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.one is not None:
        if args.one not in {r["obsid"] for r in campaign.read(output / "release_bounds.json")["obsids"]}:
            raise ValueError("Unplanned supplementary ObsID")
        scan_one(args.one, output)
        return
    (output / "observations").mkdir(exist_ok=True)
    entries = prepare(args.evidence.resolve(), output)
    deadline = time.monotonic() + 300 * 60
    for entry in entries:
        oid = entry["obsid"]
        path = output / "observations" / f"{oid}.json.gz"
        if path.is_file() and campaign.read_obs(path)["status"] == "complete":
            continue
        remaining = deadline - time.monotonic()
        if remaining < 120:
            break
        try:
            subprocess.run([sys.executable, __file__, "--one", str(oid), "--output", str(output),
                            "--evidence", str(args.evidence)], check=True, timeout=remaining)
        except (subprocess.SubprocessError, OSError) as exc:
            row = campaign.read_obs(path) if path.exists() else {"obsid":oid, "products":[], "errors":[]}
            row["status"] = "error"
            row["errors"].append(f"{type(exc).__name__}: {exc}")
            campaign.save_obs(path, row)
        summarize(output, entries)
    statuses = summarize(output, entries)
    if any(r["worker_status"] != "complete" for r in statuses):
        raise SystemExit("Supplement partially completed. See summary.json and individual input receipts.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--one", type=int, help=argparse.SUPPRESS)
    main(parser.parse_args())
