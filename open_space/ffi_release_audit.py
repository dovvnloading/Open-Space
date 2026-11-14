"""Audit a provisional FFI grid summary and hash its external supporting files.

This verifies arithmetic, acquisition-plan shape and file integrity. It does not
rerun the astronomy or turn a provisional claim into a scientific discovery.
No third-party packages are required.
"""

import argparse
from datetime import date
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from zipfile import BadZipFile, ZipFile


_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _integer(value, label, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


def _number(value, label, minimum=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{label} must be finite")
    if minimum is not None and value <= minimum:
        raise ValueError(f"{label} must be > {minimum}")
    return float(value)


def validate_report(report):
    """Return derived statistics after validating internal consistency only."""
    if report.get("schema_version") != 1:
        raise ValueError("Expected schema_version=1")
    if report.get("evidence_state") != "external_unverified":
        raise ValueError("Do not silently upgrade a provisional report")
    date.fromisoformat(report["reported_on"])
    baseline = report["baseline_commit_sha"]
    if not re.fullmatch(r"[0-9a-f]{40}", baseline):
        raise ValueError("Invalid baseline commit hash")

    plan = report["field_plan"]
    ra, dec = plan["origin_icrs_deg"]
    dra, ddec = plan["step_ra_dec_deg"]
    ra = _number(ra, "origin RA")
    dec = _number(dec, "origin Dec")
    dra = _number(dra, "RA step", 0)
    ddec = _number(ddec, "Dec step", 0)
    columns = _integer(plan["columns"], "columns", 1)
    rows = _integer(plan["rows"], "rows", 1)
    _integer(plan["size_pixels"], "size_pixels", 11)
    if plan["size_pixels"] > 41 or not 0 <= ra < 360 or not -90 < dec < 90:
        raise ValueError("Invalid sky cutout geometry")
    if dec + (rows - 1) * ddec >= 90 or dra * (columns - 1) >= 360:
        raise ValueError("Grid extends outside valid coordinate range")
    sectors = plan["sectors"]
    if len(sectors) != 2 or not all(type(x) is int for x in sectors) or not 1 <= sectors[0] < sectors[1]:
        raise ValueError("Expected two increasing TESS sectors")
    centers = [[round((ra + x*dra) % 360, 9), round(dec + y*ddec, 9)]
               for y in range(rows) for x in range(columns)]

    scan = report["reported_scan"]
    first = _integer(scan["first_sector_peaks"], "first_sector_peaks")
    outcomes = scan["first_sector_dispositions"]
    needed = {"recurrent", "image_artifact_veto", "single_sector", "frequency_inconsistent"}
    if set(outcomes) != needed:
        raise ValueError("Incomplete first-sector disposition partition")
    for k in needed:
        _integer(outcomes[k], k)
    if sum(outcomes.values()) != first:
        raise ValueError("First-sector dispositions do not sum to first-sector peaks")
    recurrent = outcomes["recurrent"]
    for k in ("later_sector_repeat", "automatic_vsx_associations", "ogle_period_associations", "direct_gaia_period_matches"):
        if _integer(scan[k], k) > recurrent:
            raise ValueError(f"{k} exceeds recurrent detections")
    area = _number(scan["two_sector_usable_area_sq_deg"], "area", 0)
    _number(scan["minimum_target_offset_arcsec"], "minimum_target_offset_arcsec", 0)

    controls = report["reported_controls"]
    injections = {}
    for item in controls["injections"]:
        amplitude = _number(item["peak_semiamplitude_e_per_s"], "amplitude", 0)
        total = _integer(item["trials"], "injection trials", 1)
        found = _integer(item["recovered"], "recovered")
        if found > total or amplitude in injections:
            raise ValueError("Invalid or duplicate injection recovery")
        injections[amplitude] = {"recovered": found, "trials": total, "fraction": found/total}
    if not injections:
        raise ValueError("No injections")
    nulls = {}
    for item in controls["nulls"]:
        method = item["method"]
        if not isinstance(method, str) or not method or method in nulls:
            raise ValueError("Invalid/duplicate null method")
        n = _integer(item["field_trials"], "null field trials", 1)
        found = _integer(item["recurrent_in_area"], "null recurring peaks")
        effective = _number(item["sq_degree_trials"], "area-trials", 0)
        nulls[method] = {"field_trials": n, "recurring_in_area": found,
                         "sq_degree_trials": effective, "recurrences_per_sq_degree_trial": found/effective}
    if not nulls:
        raise ValueError("No null controls")

    candidates = report["reported_source_followup"]
    if not candidates:
        raise ValueError("Missing follow-up cases")
    names = set()
    for item in candidates:
        name = item["catalog_source"]
        if not isinstance(name, str) or not name or name in names:
            raise ValueError("Invalid/duplicate source")
        names.add(name)
        measured = _number(item["tess_period_days_approx"], "TESS period", 0)
        period = _number(item["catalog_period_days"], "catalog period", 0)
        order = _integer(item["harmonic_order"], "harmonic_order", 1)
        if abs(period/order/measured - 1) > 0.03:
            raise ValueError(f"Source harmonic differs too far from TESS period: {name}")
        for prefix in ("blind_new_sectors", "heldout_phase"):
            value = item[prefix]
            if _integer(value["success"], prefix + ".success") > _integer(value["total"], prefix + ".total", 1):
                raise ValueError("Success count exceeds eligible observations")
        if item.get("icrs_deg") is not None:
            center = item["icrs_deg"]
            if len(center) != 2 or not 0 <= _number(center[0], "source RA") < 360 or not -90 < _number(center[1], "source Dec") < 90:
                raise ValueError("Invalid source coordinates")

    filenames = set()
    for item in report["external_artifacts"]:
        name = item["filename"]
        if not isinstance(name, str) or not name or Path(name).name != name or name in (".", "..") or "\\" in name or name in filenames:
            raise ValueError("External artifact names must be distinct safe basenames")
        filenames.add(name)
        _integer(item["bytes"], f"{name}.bytes", 1)
        if "sha256" in item and not _SHA256.fullmatch(item["sha256"]):
            raise ValueError("Invalid SHA-256 in external artifact description")

    return {"grid_centers_icrs_deg": centers, "planned_fields": rows*columns,
            "first_sector_peaks": first, "recurrent_sources": recurrent,
            "usable_area_sq_deg": area, "injection_recovery": {str(k): v for k,v in sorted(injections.items())},
            "null_controls": nulls, "described_followup_cases": len(names),
            "external_artifacts_expected": len(filenames),
            "scientific_reproduction": "not_performed"}


def audit_artifacts(report, evidence_dir=None, check_zip_crc=False):
    """Never treat file presence or hashes alone as scientific reproduction."""
    results = []
    for item in report["external_artifacts"]:
        row = {"filename": item["filename"], "expected_bytes": item["bytes"],
               "status": "not_available"}
        if evidence_dir is not None:
            path = Path(evidence_dir) / item["filename"]
            if path.is_symlink():
                row["status"] = "symlink_refused"
            elif path.is_file():
                size = path.stat().st_size
                row["observed_bytes"] = size
                if size != item["bytes"]:
                    row["status"] = "size_mismatch"
                else:
                    digest = sha256()
                    with path.open("rb") as stream:
                        for chunk in iter(lambda: stream.read(1024*1024), b""):
                            digest.update(chunk)
                    row["sha256"] = digest.hexdigest()
                    if item.get("sha256") and item["sha256"] != row["sha256"]:
                        row["status"] = "hash_mismatch"
                    else:
                        row["status"] = "locally_hashed"
                        if check_zip_crc and path.suffix.lower() == ".zip":
                            try:
                                with ZipFile(path) as archive:
                                    bad = archive.testzip()
                                if bad is not None:
                                    row["status"] = "zip_crc_failure"
                                    row["first_bad_member"] = bad
                                else:
                                    row["zip_crc_checked"] = True
                            except (BadZipFile, OSError, EOFError) as exc:
                                row["status"] = "zip_invalid"
                                row["reason"] = str(exc)
        results.append(row)
    all_present = bool(results) and all(x["status"] == "locally_hashed" for x in results)
    return {"state": "all_files_locally_hashed" if all_present else "external_evidence_incomplete",
            "scope": "file_presence_size_sha256_and_optional_zip_crc_not_scientific_reproduction",
            "artifacts": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=Path("results/ffi-lmc-east-preliminary.json"))
    parser.add_argument("--evidence-dir", type=Path)
    parser.add_argument("--check-zip-crc", action="store_true")
    parser.add_argument("--require-all", action="store_true")
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    raw = options.summary.read_bytes()
    report = json.loads(raw)
    result = {"schema_version": 1, "summary_sha256": sha256(raw).hexdigest(),
              "evidence_state": report["evidence_state"],
              "internal_consistency": validate_report(report),
              "supporting_files": audit_artifacts(report, options.evidence_dir, options.check_zip_crc)}
    content = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if options.output:
        options.output.parent.mkdir(parents=True, exist_ok=True)
        options.output.write_text(content, encoding="utf-8")
    else:
        print(content, end="")
    if options.require_all and result["supporting_files"]["state"] != "all_files_locally_hashed":
        parser.exit(2, "External evidence is incomplete; summary remains provisional.\n")


if __name__ == "__main__":
    main()
