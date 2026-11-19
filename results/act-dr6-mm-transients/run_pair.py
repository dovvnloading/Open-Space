#!/usr/bin/env python3
"""Single staged FITS-pair white-noise scan; explicitly NOT a survey detection.

A usable scientific run additionally needs measured ACT beam/scan-noise,
independent jackknife, injection completeness, multiple-trials correction,
time/position checks and a frozen full-observation denominator.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from mf_core import matched_filter, negative_tail_count, peaks
from freeze_mapcat import digest


def read_image(path):
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(path, memmap=True) as hdus:
        head = hdus[0].header
        wcs = WCS(head).celestial
        data = np.asarray(hdus[0].data, dtype=np.float64)
        if data.ndim == 3 and data.shape[0] == 3:
            data = data[0]  # T of T,Q,U maps
        if data.ndim != 2:
            raise ValueError(f"Unexpected FITS image dimensions: {path}")
        return data.copy(), wcs.to_header(relax=True).tostring(sep="\\n", padding=False), wcs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("map", "ivar", "reference", "reference-ivar", "psf"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--obs-id", required=True, help="Exact frozen MapCat map_id")
    p.add_argument("--reference-receipt", type=Path, required=True,
                   help="JSON metadata proving leave-one-out exclusion of obs-id")
    p.add_argument("--beam-receipt", type=Path, required=True,
                   help="JSON with measured ACT beam provenance and psf_sha256")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--threshold", type=float, default=7.0)
    args = p.parse_args()
    ref_receipt = json.loads(args.reference_receipt.read_text())
    if args.obs_id not in ref_receipt.get("excluded_map_ids", []):
        p.error("Reference has not documented exclusion of observation")
    beam = json.loads(args.beam_receipt.read_text())
    if beam.get("kind") != "measured_ACT_beam" or beam.get("psf_sha256") != digest(args.psf):
        p.error("Missing or mismatched measured ACT beam receipt")
    psf = np.load(args.psf, allow_pickle=False)
    data = [read_image(getattr(args, attr.replace("-", "_"))) for attr in
            ("map", "reference", "ivar", "reference-ivar")]
    grid_signature = data[0][1]
    if any(a[1] != grid_signature for a in data):
        p.error("FITS WCS does not match exactly: reproject with recorded interpolation first")
    if any(a[0].shape != data[0][0].shape for a in data):
        p.error("FITS shapes differ")
    obs, reference, ivar, refivar = [a[0] for a in data]
    result = matched_filter(obs, reference, ivar, refivar, psf, independent_reference=True)
    positions = peaks(result["snr_white"], args.threshold)
    celestial_wcs = data[0][2]
    for found in positions:
        found["ra_deg"], found["dec_deg"] = [
            float(v) for v in celestial_wcs.all_pix2world(found["x"], found["y"], 0)]
        found["amp_native"] = float(result["amplitude_native"][found["y"], found["x"]])
        found["sigma_white_native"] = float(result["sigma_white_native"][found["y"], found["x"]])
        found["disposition"] = "uncalibrated_pixel_peak_not_a_transient"
    payload = {
        "status": "uncalibrated_white_noise_screen_only",
        "utc": datetime.now(timezone.utc).isoformat(),
        "obs_id": args.obs_id,
        "files_sha256": {name: digest(getattr(args, name.replace("-", "_")))
                         for name in ("map", "ivar", "reference", "reference-ivar", "psf")},
        "reference_receipt_sha256": digest(args.reference_receipt),
        "beam_receipt_sha256": digest(args.beam_receipt),
        "threshold_white_snr": args.threshold,
        "valid_centres": int(result["eligible"].sum()),
        "white_positive_peaks": positions,
        "white_negative_pixels": negative_tail_count(result["snr_white"], args.threshold),
        "noise_correlations_calibrated": False,
        "jackknife_completed": False,
        "source_injection_completeness_completed": False,
        "catalog_and_ephemeris_vetting_completed": False,
        "trial_factor_calibrated": False,
        "detection_claim": False,
        "valid_rate_limit": False,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\\n")
    print(json.dumps({k: v for k, v in payload.items() if k != "white_positive_peaks"}, indent=2))
    print("Uncalibrated peaks:", len(positions))


if __name__ == "__main__":
    main()
