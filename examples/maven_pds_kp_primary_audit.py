#!/usr/bin/env python3
"""Audit the 2025-08-24 Mars energetic-particle peak against ORIGINAL NASA PDS KP.

Independently obtain the NASA PDS4 2025-08-24 v26/r01 TAB and XML label, verify
the label's official TAB MD5 digest, decode field widths/offsets from that PDS4
label (NOT the Parquet schema), and compare raw spacecraft measurements with
the separate third-party LASP-KP Parquet mirror used in the full survey.

An unsuccessful network fetch or any failed numerical assertion is NOT a
validated result. No cross-product independent calibration is claimed: these
are two representations of MAVEN KP; dedicated SEP L2 remains independent.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import tempfile
import time
import urllib.error
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

DATASET = "maven.insitu.calibrated/data/2025/08/mvn_kp_insitu_20250824_v26_r01"
OFFICIAL_MD5 = "377f757a180b34c7af46b2937831c181"
OFFICIAL_BYTES = 45698463
SOURCE = "NASA PDS Planetary Plasma Interactions, KP collection 41.0, version v26_r01 (2026-06-10)"
HF_SOURCE = "https://huggingface.co/datasets/juliensimon/nasa-maven-kp-insitu/resolve/main/data/sep/sep.parquet"
SIDES = ("1-F", "1-R", "2-F", "2-R")
CHANNELS = {
    **{
        f"sep_ion_fov{x.lower().replace('-', '')}": f"SEP:Ion Flux (30-1000 keV), FOV {x}"
        for x in SIDES
    },
    **{
        f"sep_electron_fov{x.lower().replace('-', '')}": f"SEP:Electron Flux (30 keV - 300 keV) - FOV {x}"
        for x in SIDES
    },
}
CHANNELS.update({k + "_quality": v + " uncertainty" for k, v in list(CHANNELS.items())})


def pds_url(extension: str) -> str:
    identifier = "pds://PPI/" + DATASET + "." + extension
    return "https://search-pdsppi.igpp.ucla.edu/ditdos/download?id=" + urllib.parse.quote(identifier, safe="")


def download(url: str, path: Path) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "Open-Space-scientific-validation/1.0"})
    for attempt in range(3):
        dig = hashlib.sha256()
        md5 = hashlib.md5()
        try:
            with urllib.request.urlopen(req, timeout=240) as response, path.open("wb") as output:
                for chunk in iter(lambda: response.read(1024 * 1024), b""):
                    dig.update(chunk)
                    md5.update(chunk)
                    output.write(chunk)
            return {"url": url, "bytes": path.stat().st_size, "sha256": dig.hexdigest(), "md5": md5.hexdigest()}
        except urllib.error.HTTPError as exc:
            if exc.code not in (429, 500, 502, 503, 504) or attempt == 2:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def item(parent, name: str) -> str:
    element = parent.find("{*}" + name)
    if element is None or element.text is None:
        raise ValueError("Missing primary-PDS label field " + name)
    return element.text.strip()


def parse_pds(path: Path, label: Path) -> tuple[pd.DataFrame, dict]:
    root = ET.parse(label).getroot()
    table = root.find(".//{*}Table_Character")
    if table is None:
        raise ValueError("PDS XML contains no Table_Character")
    start, records = int(item(table, "offset")), int(item(table, "records"))
    record = table.find("{*}Record_Character")
    if record is None:
        raise ValueError("PDS XML missing Record_Character")
    length = int(item(record, "record_length"))
    all_fields = {}
    for field in record.findall("{*}Field_Character"):
        all_fields[item(field, "name")] = {
            "start": int(item(field, "field_location")) - 1,
            "length": int(item(field, "field_length")),
        }
    required = {"time": "Time (UTC/SCET)", **CHANNELS}
    missing = sorted(set(required.values()) - set(all_fields))
    if missing:
        raise ValueError("Missing expected primary KP SEP field names: " + str(missing))
    raw = path.read_bytes()
    if len(raw) != start + records * length:
        raise ValueError(f"PDS label record length inconsistent: {len(raw)} vs {start + records*length}")
    selected = {alias: all_fields[label_name] for alias, label_name in required.items()}
    rows = []
    for idx in range(records):
        row = raw[start + idx * length : start + (idx + 1) * length]
        rec = {}
        for alias, spec in selected.items():
            val = row[spec["start"] : spec["start"] + spec["length"]].decode("ascii").strip()
            if alias == "time":
                rec[alias] = val
            else:
                try:
                    rec[alias] = float(val)
                except ValueError:
                    rec[alias] = np.nan
        rows.append(rec)
    data = pd.DataFrame(rows)
    data["time"] = pd.to_datetime(data["time"], errors="coerce", utc=True)
    if data["time"].isna().any():
        raise ValueError("PDS primary source contains unparseable times")
    return data, {
        "expected_records": records, "parsed_records": int(len(data)),
        "record_length": length, "data_offset": start,
        "label_columns_used": required,
    }


def audit(output: Path) -> dict:
    report = {"status": "not yet validated", "source": SOURCE,
              "official_primary_pds_id": DATASET,
              "official_primary_md5": OFFICIAL_MD5,
              "source_link": "https://search-pdsppi.igpp.ucla.edu/search/view/?f=no&id=pds%3A%2F%2FPPI%2F" + DATASET.replace("/", "%2F") + "&o=1",
              "limitations": ["Different file conversions of the SAME underlying KP sensor data; not an independent instrument or an independent SEP L2 calibration.", "This checks the high-energy SEP peak day, not full-survey completeness or all named CMEs."]}
    try:
        with tempfile.TemporaryDirectory(prefix="maven-pds-audit-") as tmp:
            root = Path(tmp)
            tab = root / "primary.tab"
            xml = root / "primary.xml"
            parquet = root / "sep.parquet"
            for ext, dest in (("xml", xml), ("tab", tab)):
                report["primary_" + ext] = download(pds_url(ext), dest)
            if report["primary_tab"]["md5"] != OFFICIAL_MD5 or report["primary_tab"]["bytes"] != OFFICIAL_BYTES:
                raise AssertionError("Primary PDS .tab failed official label MD5/size check")
            pds, schema = parse_pds(tab, xml)
            report["pds_schema"] = schema
            report["mirror_parquet"] = download(HF_SOURCE, parquet)
            if report["mirror_parquet"]["sha256"] != "e3ea044212db893dcdcc050c5b1aa28257b05f483c9952258bc715e9631dec03":
                raise AssertionError("Parquet mirror changed since locked full-survey input")
            mirror = pd.read_parquet(parquet, columns=["time"] + list(CHANNELS))
            mirror["time"] = pd.to_datetime(mirror["time"], utc=True)
            date_start, date_end = pd.Timestamp("2025-08-24", tz="UTC"), pd.Timestamp("2025-08-25", tz="UTC")
            mirror = mirror[(mirror["time"] >= date_start) & (mirror["time"] < date_end)]
            peak_start = pd.Timestamp("2025-08-24T05:00:00Z")
            peak_end = peak_start + pd.Timedelta(hours=1)
            mirror_peak = mirror[(mirror["time"] >= peak_start) & (mirror["time"] < peak_end)]
            pds_peak = pds[(pds["time"] >= peak_start) & (pds["time"] < peak_end)]
            both = mirror_peak.merge(pds_peak, on="time", how="inner", suffixes=("_mirror", "_pds"), validate="one_to_one")
            report["overlap"] = {
                "pds_full_day_rows": int(len(pds)),
                "mirror_full_day_rows": int(len(mirror)),
                "mirror_peak_hour_rows": int(len(mirror_peak)),
                "pds_peak_hour_rows": int(len(pds_peak)),
                "matched_peak_hour_times": int(len(both)),
                "fraction_mirror_peak_rows_matched": float(len(both) / len(mirror_peak)) if len(mirror_peak) else 0.0,
            }
            if len(both) < 300 or len(both) / max(1, len(mirror_peak)) < 0.90:
                raise AssertionError("Insufficient 1-to-1 time matching to primary PDS")
            deltas = {}
            for column in CHANNELS:
                a = both[column + "_mirror"].to_numpy(dtype=float)
                b = both[column + "_pds"].to_numpy(dtype=float)
                valid = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0) & (a < 1e10) & (b < 1e10)
                if valid.sum() < 200:
                    raise AssertionError(f"Too few comparable high-flux PDS samples for {column}")
                diff = np.abs(a[valid] - b[valid]) / np.maximum(1, np.abs(b[valid]))
                deltas[column] = {
                    "matched_positive_values": int(valid.sum()),
                    "median_absolute_fractional_difference": float(np.median(diff)),
                    "p95_absolute_fractional_difference": float(np.quantile(diff, 0.95)),
                }
            report["flux_and_uncertainty_checks"] = deltas
            # Comparability of independent archive encoding is strict; if the PDS
            # v26 and the mirror differ, report that fact, never silently promote.
            if max(v["median_absolute_fractional_difference"] for v in deltas.values()) > 0.05:
                raise AssertionError("Primary NASA PDS and mirror flux differ by >5% median in one or more channels")
            report["status"] = "passed: original NASA PDS Aug24 KP flux and 1-sigma uncertainty values independently retrieved and numerically consistent with mirror"
    except Exception as exc:
        report["status"] = "NOT validated: " + str(exc)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + "\n")
        raise
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=Path("results/maven-2025-sep-events/primary_pds_audit.json"))
    args = p.parse_args()
    result = audit(args.output)
    print(json.dumps({"status": result["status"], "overlap": result["overlap"], "worst_median_fractional_diff": max(x["median_absolute_fractional_difference"] for x in result["flux_and_uncertainty_checks"].values())}, indent=2))
