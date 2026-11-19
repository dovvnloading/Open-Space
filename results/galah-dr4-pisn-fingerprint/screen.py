#!/usr/bin/env python3
"""GALAH DR4 full-file, flag-aware preregistered PISN chemical screen.

Requires numpy and astropy. No candidate is a supernova model match.
The published official 240705 allstar FITS is downloaded unmodified; the
complete input SHA-256 and byte count are recorded. Never use flag_fe_h.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
import urllib.request
from datetime import datetime, timezone
import numpy as np

HERE = Path(__file__).resolve().parent
ELTS = ("na", "al", "mg", "co", "zn", "ba", "sr")
FIVE = ELTS[:5]
URL = "https://cloud.datacentral.org.au/teamdata/GALAH/public/GALAH_DR4/catalogs/galah_dr4_allstar_240705.fits"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def put_json(path: Path, data) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def download(path: Path) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.stat().st_size > 100_000_000:
        return {"source": "existing_local_cache", "url": URL, "bytes": path.stat().st_size, "sha256": sha256(path)}
    last_error = None
    for attempt in range(1, 4):
        part = path.with_name(path.name + ".part")
        part.unlink(missing_ok=True)
        try:
            req = urllib.request.Request(URL, headers={"User-Agent": "Open-Space-GALAH-DR4-Research/1.0"})
            with urllib.request.urlopen(req, timeout=150) as resp, part.open("wb") as out:
                headers = {
                    "content_length": resp.headers.get("Content-Length"),
                    "last_modified": resp.headers.get("Last-Modified"),
                    "etag": resp.headers.get("ETag"),
                    "response_url": resp.geturl(),
                    "status": resp.status,
                }
                h = hashlib.sha256()
                nbytes = 0
                while block := resp.read(8 * 1024 * 1024):
                    out.write(block)
                    h.update(block)
                    nbytes += len(block)
                if nbytes < 100_000_000:
                    raise ValueError(f"Unusually small FITS download: {nbytes} bytes")
                part.replace(path)
                return {"source": "live_download", "url": URL, "attempt": attempt,
                        "bytes": nbytes, "sha256": h.hexdigest(), "http": headers}
        except Exception as exc:
            last_error = repr(exc)
            print(f"Download attempt {attempt}/3 failed: {last_error}", flush=True)
            part.unlink(missing_ok=True)
            if attempt < 3:
                time.sleep(4 * attempt)
    raise RuntimeError(f"Official FITS unavailable after three attempts: {last_error}")


def gate_masks(abund: dict[str, np.ndarray], cuts: dict) -> tuple[np.ndarray, np.ndarray]:
    """Thresholds on *measured* values; any unmeasured value must be masked outside."""
    five = ((abund["na"] - abund["mg"] <= cuts["na_minus_mg_max"])
            & (abund["al"] - abund["mg"] <= cuts["al_minus_mg_max"])
            & (abund["co"] <= cuts["co_fe_max"])
            & (abund["zn"] <= cuts["zn_fe_max"]))
    seven = five & (abund["ba"] <= cuts["ba_fe_max"]) & (abund["sr"] <= cuts["sr_fe_max"])
    return five, seven


def injection(mask: np.ndarray, err: dict[str, np.ndarray], criteria: dict,
              include_neutron: bool, rng: np.random.Generator) -> dict:
    idx = np.flatnonzero(mask)
    trials = int(criteria["injection_trials"])
    if not len(idx):
        return {"eligible_stars": 0, "trials_per_star": trials, "passing_trials": 0,
                "completeness_conditional_on_flags": None, "uncertainty": None}
    truth = criteria["injection_dex"]
    passed = 0
    cuts = criteria["gate"]
    for _ in range(trials):
        draws = {k: truth[k] + rng.normal(size=len(idx)) * err[k][idx]
                 for k in (ELTS if include_neutron else FIVE)}
        matched = ((draws["na"] - draws["mg"] <= cuts["na_minus_mg_max"])
                   & (draws["al"] - draws["mg"] <= cuts["al_minus_mg_max"])
                   & (draws["co"] <= cuts["co_fe_max"])
                   & (draws["zn"] <= cuts["zn_fe_max"]))
        if include_neutron:
            matched &= (draws["ba"] <= cuts["ba_fe_max"]) & (draws["sr"] <= cuts["sr_fe_max"])
        passed += int(matched.sum())
    n = len(idx) * trials
    p = passed / n
    return {"eligible_stars": int(len(idx)), "trials_per_star": trials,
            "passing_trials": passed, "completeness_conditional_on_flags": p,
            "monte_carlo_standard_error": math.sqrt(p * (1-p)/n),
            "warning": "Measured-abundance flags are not independent of true composition; this is not astrophysical search completeness."}


def scan(table, criteria: dict, input_receipt: dict) -> tuple[dict, list[dict], list[dict]]:
    names = set(table.columns.names)
    required = {"sobject_id", "gaiadr3_source_id", "ra", "dec",
                "fe_h", "e_fe_h", "flag_sp", "snr_px_ccd3"}
    for el in ELTS:
        required |= {f"{el}_fe", f"e_{el}_fe", f"flag_{el}_fe"}
    missing = required - names
    if missing:
        raise ValueError(f"Missing columns in official allstar: {sorted(missing)}")
    def num(name):
        return np.asarray(table[name], dtype=np.float64)
    def integer(name):
        return np.asarray(table[name], dtype=np.int64)
    n = len(table)
    if n != 917588:
        raise ValueError(f"Expected 917588 rows in frozen release; got {n}")
    feh, efeh, snr, fsp = num("fe_h"), num("e_fe_h"), num("snr_px_ccd3"), integer("flag_sp")
    zcut = criteria["metallicity_dex"]
    qcut = criteria["quality"]
    metal = np.isfinite(feh) & (feh >= zcut["min"]) & (feh <= zcut["max"])
    quality = ((fsp == qcut["flag_sp"]) & np.isfinite(snr) & (snr > qcut["snr_px_ccd3_gt"])
               & np.isfinite(efeh) & (efeh > 0))
    base = metal & quality
    avail, abund, err, qc_counts = {}, {}, {}, {}
    for el in ELTS:
        v, e, flag = num(f"{el}_fe"), num(f"e_{el}_fe"), integer(f"flag_{el}_fe")
        a = ((flag == qcut["per_element_flag"]) & np.isfinite(v)
             & np.isfinite(e) & (e > qcut["per_element_sigma_dex_min_exclusive"])
             & (e <= qcut["per_element_sigma_dex_max"]))
        abund[el], err[el], avail[el] = v, e, a
        qc_counts[el] = {
            "unflagged_finite_uncertainty_eligible_in_base": int((base & a).sum()),
            "flag_zero_in_base": int((base & (flag == 0)).sum()),
            "flag_upper_limit_bit_in_base": int((base & ((flag & 1) != 0)).sum()),
            "flag_undetectable_bit_in_base": int((base & ((flag & 2) != 0)).sum()),
            "missing_or_rejected_in_base": int((base & ~a).sum())
        }
    mask5 = base.copy()
    for el in FIVE:
        mask5 &= avail[el]
    mask7 = mask5 & avail["ba"] & avail["sr"]
    raw5, raw7 = gate_masks(abund, criteria["gate"])
    hit5 = mask5 & raw5
    hit7 = mask7 & raw7
    odd = base & avail["na"] & avail["al"] & avail["mg"]
    iron = base & avail["co"] & avail["zn"]
    threshold_odd = odd & ((abund["na"]-abund["mg"]) <= criteria["gate"]["na_minus_mg_max"]) & ((abund["al"]-abund["mg"]) <= criteria["gate"]["al_minus_mg_max"])
    rng = np.random.default_rng(int(criteria["seed"]))
    pow5 = injection(mask5, err, criteria, False, rng)
    pow7 = injection(mask7, err, criteria, True, rng)
    ids = integer("sobject_id")
    gaia = integer("gaiadr3_source_id")
    ra, dec = num("ra"), num("dec")
    def row(i):
        output = {"sobject_id": int(ids[i]), "gaia_dr3_source_id": int(gaia[i]),
                  "ra_deg": float(ra[i]), "dec_deg": float(dec[i]),
                  "fe_h": float(feh[i]), "e_fe_h": float(efeh[i]),
                  "snr_px_ccd3": float(snr[i]), "five_gate": bool(hit5[i]),
                  "seven_gate": bool(hit7[i])}
        for el in ELTS:
            output[f"{el}_fe"] = float(abund[el][i]) if avail[el][i] else None
            output[f"e_{el}_fe"] = float(err[el][i]) if avail[el][i] else None
        if avail["na"][i] and avail["mg"][i]:
            output["na_minus_mg"] = float(abund["na"][i]-abund["mg"][i])
        else:
            output["na_minus_mg"] = None
        if avail["al"][i] and avail["mg"][i]:
            output["al_minus_mg"] = float(abund["al"][i]-abund["mg"][i])
        else:
            output["al_minus_mg"] = None
        return output
    # Preserve every complete five-element diagnostic, even when it fails the gate.
    complete5_rows = [row(int(i)) for i in np.flatnonzero(mask5)]
    hits = [row(int(i)) for i in np.flatnonzero(hit5)]
    summary = {
        "data_receipt": input_receipt,
        "source_model":"GALAH DR4 official 240705 allstar; one row per star",
        "rows": {"catalogue": n, "fe_h_range_only": int(metal.sum()),
                 "overall_flag_and_SNR_only": int(quality.sum()),
                 "metal_poor_quality_base": int(base.sum()),
                 "odd_even_3_measured": int(odd.sum()),
                 "odd_even_3_pass_thresholds": int(threshold_odd.sum()),
                 "co_and_zn_measured": int(iron.sum()),
                 "all_five_measured": int(mask5.sum()),
                 "all_seven_measured": int(mask7.sum()),
                 "five_pattern_pass": int(hit5.sum()),
                 "seven_pattern_pass": int(hit7.sum())},
        "element_coverage": qc_counts,
        "injection":{"assumed_true_dex": criteria["injection_dex"],
                     "five": pow5, "seven": pow7,
                     "five_total_base_completeness_proxy": (len(complete5_rows) / max(1,int(base.sum())) * pow5["completeness_conditional_on_flags"]) if len(complete5_rows) else 0.0,
                     "seven_total_base_completeness_proxy": (int(mask7.sum())/max(1,int(base.sum()))*pow7["completeness_conditional_on_flags"]) if int(mask7.sum()) else 0.0},
        "interpretation": "These hard-threshold flags do NOT establish PISN origin. Published PISN/CCSN/Ia yields, censored-data likelihoods, survey cross-checks and line vetting remain uncompleted.",
        "criteria_sha256": sha256(HERE / "criteria.json"),
        "script_sha256": sha256(Path(__file__)),
        "run_utc": datetime.now(timezone.utc).isoformat(),
    }
    return summary, complete5_rows, hits


def export_csv(path: Path, rows: list[dict]) -> None:
    cols = ("sobject_id", "gaia_dr3_source_id", "ra_deg", "dec_deg", "fe_h", "e_fe_h",
            "snr_px_ccd3", "five_gate", "seven_gate", "na_minus_mg", "al_minus_mg")
    cols += tuple(f"{el}_fe" for el in ELTS)
    cols += tuple(f"e_{el}_fe" for el in ELTS)
    with path.open("w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=cols)
        wr.writeheader()
        wr.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=None)
    ap.add_argument("--cache", type=Path, default=Path("galah_dr4_allstar_240705.fits"))
    ap.add_argument("--out", type=Path, default=HERE / "run")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    try:
        crit = json.loads((HERE / "criteria.json").read_text())
        receipt = download(args.cache) if args.input is None else {
            "source": "user_supplied_file", "path": str(args.input), "url": URL,
            "bytes": args.input.stat().st_size, "sha256": sha256(args.input)}
        path = args.cache if args.input is None else args.input
        from astropy.io import fits
        with fits.open(path, memmap=True, lazy_load_hdus=True) as hdus:
            data = hdus[1].data
            if data is None:
                raise ValueError("No FITS table at HDU 1")
            summary, full, hits = scan(data, crit, receipt)
        put_json(args.out / "scan_summary.json", summary)
        export_csv(args.out / "complete_five_rows.csv", full)
        export_csv(args.out / "five_gate_survivors.csv", hits)
        print("GALAH_RESULT", json.dumps(summary["rows"], sort_keys=True), flush=True)
        print("GALAH_INPUT_SHA256", receipt["sha256"], flush=True)
    except Exception as exc:
        put_json(args.out / "scan_blocked.json", {"status":"blocked", "error":repr(exc),
                 "utc": datetime.now(timezone.utc).isoformat(),
                 "script_sha256":sha256(Path(__file__))})
        print(f"GALAH_SCAN_BLOCKED {exc!r}", file=sys.stderr, flush=True)
        raise


if __name__ == "__main__":
    main()
