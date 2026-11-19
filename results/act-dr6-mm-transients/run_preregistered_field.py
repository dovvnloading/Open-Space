#!/usr/bin/env python3
"""Bounded, preregistered ACT pa4/f150 real-FITS ROI attempt.

All failures are retained. Remote FITS ranges are not full-product checksums.
Even successful outputs are a non-discovery, non-rate, beam-transfer-unvalidated
pilot; no Horizons/catalog verdict can be inferred from its SNR peaks.
"""
from __future__ import annotations
import gzip
import hashlib
import io
import json
import math
import sys
import tarfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from scipy.special import j0

from freeze_mapcat import ROOT, digest
from mf_core import matched_filter, peaks, inject, gaussian_peak_psf

BASE = Path(__file__).resolve().parent
MANIFEST = BASE / "eligible-act-depth1.jsonl.gz"
MANIFEST_SHA = "b6a0c020e3b03bc23aac66536ae1b31cb8889cbc3295c1d2acf16f8de3952503"
DATE_START, DATE_END = "2018-09-10T00:00:00", "2018-09-13T00:00:00"
REF_START, REF_END = "2018-08-21T00:00:00", "2018-10-04T00:00:00"
CENTER_RA, CENTER_DEC, RADIUS_DEG = 301.9952, 16.1652, 0.30
PIX_DEG, N = 0.5 / 60, 161
MAX_HEADER = 65536
MAX_WINDOW = 500
MAX_REQUESTS_PER_FITS = 520
USER_AGENT = "Open-Space-ACT-fixed-f150-science-pilot/1.0"
# Published f150 AdvACT resolution, not a measured pa4 depth-1 transfer function.
# +/-0.2 arcmin is an intentionally broad ASSUMED sensitivity bracket, NOT
# a measured 1-sigma error. Every screen and injection tests all three widths.
BEAM_FWHM_ARCMIN = (1.2, 1.4, 1.6)
BEAM_REFERENCE = "https://doi.org/10.3847/1538-4357/ace599"
BEAM_ARCHIVE_NOT_APPLICABLE = "https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/act_dr6.02_harmonic_beams_profiles_get.html"
CACHE_FILE = BASE / "roi-footprint-cache.json"
KNOWN_RA, KNOWN_DEC = 301.9952, 16.1652



def in_interval(value, low, high):
    # MapCat start_time uses a SPACE rather than the ISO 'T' separator.
    # String comparisons against a 'T' boundary silently excluded Sept 10,
    # incorrectly allowed Sept 10 as a reference (the original run).
    return datetime.fromisoformat(low) <= datetime.fromisoformat(str(value)) < datetime.fromisoformat(high)


def selection():
    if digest(MANIFEST) != MANIFEST_SHA:
        raise ValueError("Frozen compressed catalog hash mismatch")
    rows = [json.loads(line) for line in gzip.open(MANIFEST, "rt", encoding="utf-8")]
    rows = [r for r in rows if r["tube_slot"] == "pa4" and r["frequency"] == "f150"]
    tests = [r for r in rows if in_interval(r["start_time"], DATE_START, DATE_END)]
    refs = [r for r in rows if in_interval(r["start_time"], REF_START, REF_END)
            and not in_interval(r["start_time"], DATE_START, DATE_END)]
    if len({r["map_id"] for r in tests + refs}) != len(tests + refs):
        raise ValueError("Duplicate map id in frozen subset")
    assert all(not in_interval(r["start_time"], DATE_START, DATE_END) for r in refs)
    return tests, refs


def load_footprint_cache():
    data = json.loads(CACHE_FILE.read_text())
    if data["frozen_manifest_sha256"] != MANIFEST_SHA or data["baseline_commit_sha"] != "d66894f7e0ad434d7d5b73d12df95a18d8fcf9f4":
        raise ValueError("Cached header-footprint provenance does not match frozen catalog/baseline")
    return data["outside_header_records"]


def digest_bytes(data):
    return hashlib.sha256(data).hexdigest()


def fetch_range(url, start, stop, *, header=False):
    if start < 0 or stop < start or stop - start + 1 > 2**20:
        raise ValueError("Invalid bounded range")
    req = Request(url, headers={"Range": f"bytes={start}-{stop}", "User-Agent": USER_AGENT})
    err = None
    for trial in range(3):
        try:
            with urlopen(req, timeout=28) as r:
                status = int(r.status)
                data = r.read(stop - start + 2)
                if status != 206:
                    if header and status == 200:
                        # Read only bounded header even if range ignored.
                        return data[:stop-start+1], {"status": status, "range": None}
                    raise ValueError(f"HTTP {status}: byte range not honored")
                expected = f"bytes {start}-{stop}/"
                cr = r.headers.get("Content-Range", "")
                if not cr.startswith(expected) or len(data) != stop-start+1:
                    raise ValueError(f"Incorrect Content-Range or length: {cr} {len(data)}")
                return data, {"status": status, "range": cr}
        except Exception as e:
            err = e
            if trial < 2:
                time.sleep(0.3 * (trial+1))
    raise RuntimeError(f"Range request failed: {type(err).__name__}: {err}")


def fits_header(url):
    raw, meta = fetch_range(url, 0, MAX_HEADER-1, header=True)
    if not raw.startswith(b"SIMPLE  ="):
        raise ValueError("Not FITS SIMPLE primary")
    cards = [raw[i:i+80] for i in range(0, len(raw)//80*80, 80)]
    end = next((i for i, card in enumerate(cards) if card[:8] == b"END     "), None)
    if end is None:
        raise ValueError("FITS END not in first 64 KiB")
    head_len = ((end + 1) * 80 + 2879)//2880*2880
    h = fits.Header.fromstring(b"".join(cards[:end+1]).decode("ascii"), sep="")
    if int(h.get("NAXIS", 0)) not in (2, 3) or int(h.get("BITPIX", 0)) not in (-32, -64):
        raise ValueError("FITS primary is not expected floating image")
    info = {
        "header_64k_sha256": digest_bytes(raw),
        "status": meta["status"],
        "content_range": meta["range"],
        "naxis": h["NAXIS"], "shape": [int(h["NAXIS2"]), int(h["NAXIS1"])],
        "bitpix": h["BITPIX"],
        "header_bytes": head_len,
    }
    return h, info


def tile_wcs():
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [CENTER_RA, CENTER_DEC]
    w.wcs.crpix = [(N+1)/2, (N+1)/2]
    w.wcs.cdelt = [-PIX_DEG, PIX_DEG]
    yy, xx = np.mgrid[:N, :N]
    ra, dec = w.all_pix2world(xx, yy, 0)
    # Exact great-circle angular aperture, including the 0.5 arcmin grid.
    sep = np.degrees(np.arccos(np.clip(
        np.sin(np.radians(dec))*math.sin(math.radians(CENTER_DEC)) +
        np.cos(np.radians(dec))*math.cos(math.radians(CENTER_DEC))*
        np.cos(np.radians(ra-CENTER_RA)), -1, 1)))
    return ra, dec, sep <= RADIUS_DEG


def access_tile(path, ra, dec, *, ivar=False):
    if not path or "://" in path or ".." in path or "?" in path or "#" in path:
        raise ValueError("Invalid frozen product path")
    url = ROOT + "/" + path.lstrip("/")
    h, meta = fits_header(url)
    w = WCS(h).celestial
    x, y = w.all_world2pix(ra, dec, 0)
    nx, ny = int(h["NAXIS1"]), int(h["NAXIS2"])
    within = np.isfinite(x) & np.isfinite(y) & (x >= 1) & (y >= 1) & (x < nx-2) & (y < ny-2)
    if not np.any(within):
        return np.full(ra.shape, np.nan), {**meta, "roi_overlap": False, "science_bytes": 0}
    if not np.all(np.isfinite(x[within])) or not np.all(np.isfinite(y[within])):
        raise ValueError("Invalid WCS transformed ROI")
    x0 = max(0, int(np.floor(np.min(x[within]))) - 1)
    x1 = min(nx, int(np.ceil(np.max(x[within]))) + 2)
    y0 = max(0, int(np.floor(np.min(y[within]))) - 1)
    y1 = min(ny, int(np.ceil(np.max(y[within]))) + 2)
    if x1-x0 > MAX_WINDOW or y1-y0 > MAX_WINDOW:
        raise ValueError("Unexpected ROI WCS pixel extent")
    if (y1-y0) > MAX_REQUESTS_PER_FITS:
        raise ValueError("ROI exceeds range-request bound")
    bpp = abs(int(h["BITPIX"]))//8
    dtype = ">f4" if bpp == 4 else ">f8"
    tile = np.empty((y1-y0, x1-x0), dtype=np.float64)
    def row(yi):
        start = int(meta["header_bytes"]) + (yi*nx+x0)*bpp
        raw, _ = fetch_range(url, start, start+(x1-x0)*bpp-1)
        return yi, raw
    records = []
    with ThreadPoolExecutor(max_workers=10) as ex:
        for yi, raw in ex.map(row, range(y0, y1)):
            tile[yi-y0] = np.frombuffer(raw, dtype=dtype).astype(np.float64)
            records.append((yi, digest_bytes(raw)))
    xs = x-x0
    ys = y-y0
    safe_x = np.where(within, xs, -100)
    safe_y = np.where(within, ys, -100)
    if ivar:
        # Propagate WHITE diagonal variances through bilinear reprojection,
        # rather than interpolating inverse variance as if it were a scalar.
        # Neighbor-to-neighbor scan covariance remains unmodelled.
        ix = np.clip(np.floor(safe_x).astype(int), 0, tile.shape[1]-2)
        iy = np.clip(np.floor(safe_y).astype(int), 0, tile.shape[0]-2)
        fx = np.clip(safe_x-ix, 0., 1.)
        fy = np.clip(safe_y-iy, 0., 1.)
        var = np.zeros(ra.shape, dtype=float)
        ok = within.copy()
        for dy in (0, 1):
            for dx in (0, 1):
                q = tile[iy+dy, ix+dx]
                a = (fy if dy else 1-fy)*(fx if dx else 1-fx)
                ok &= (q > 0) & np.isfinite(q)
                var += np.divide(a*a, q, out=np.full(ra.shape, np.inf),
                                 where=(q>0)&np.isfinite(q))
        interp = np.divide(1., var, out=np.zeros_like(var), where=ok & np.isfinite(var) & (var>0))
        interp[~ok] = 0
    else:
        interp = map_coordinates(tile, [safe_y, safe_x], order=1, mode="constant", cval=np.nan)
        interp[~within] = np.nan
    meta.update({
        "roi_overlap": True, "url": url,
        "science_bytes": (x1-x0)*(y1-y0)*bpp,
        "range_requests": len(records), "window": [x0, x1, y0, y1],
        "ordered_range_digests_sha256": digest_bytes(
            "".join(sha for _, sha in sorted(records)).encode()),
        "finite_roi_pixels": int(np.isfinite(interp).sum()),
    })
    return interp, meta


def read_pair(row, ra, dec, aperture, cache=None):
    rec = {"map_id": str(row["map_id"]), "map_name": row["map_name"],
           "start_time": row["start_time"], "stop_time": row["stop_time"]}
    cached = (cache or {}).get(str(row["map_id"]))
    if cached:
        if (cached["map_name"] != row["map_name"] or cached["start_time"] != row["start_time"] or
                cached["stop_time"] != row["stop_time"] or cached["status"] != "outside_fits_footprint"):
            raise ValueError("Frozen row disagrees with baseline science-run FITS-header exclusion")
        # Reuse ONLY previously WCS-verified noncoverage, never science pixels.
        # An origin product mutation would invalidate this optimization.
        return {**cached, "reuse": "prior_run_FITS_header_outside_only_not_rechecked"}, None
    try:
        m, rm = access_tile(row["map_path"], ra, dec)
        i, ri = access_tile(row["ivar_path"], ra, dec, ivar=True)
        if not rm.get("roi_overlap") or not ri.get("roi_overlap"):
            rec["status"] = "outside_fits_footprint"
            rec["map_header"] = rm
            rec["ivar_header"] = ri
            return rec, None
        i = np.where(np.isfinite(i) & (i > 0), i, 0)
        valid = aperture & np.isfinite(m) & (i > 0)
        rec.update({"status": "science_pixels_read", "map_header": rm, "ivar_header": ri,
                    "valid_aperture_pixels": int(valid.sum()), "native_map_unit": "uK_CMB_expected_unconfirmed",
                    "map_median_valid_native": float(np.median(m[valid])) if np.any(valid) else None,
                    "ivar_median_valid_native": float(np.median(i[valid])) if np.any(valid) else None})
        if not np.any(valid):
            rec["status"] = "no_positive_ivar_in_aperture"
            return rec, None
        return rec, (np.where(valid, m, 0.), np.where(valid, i, 0.))
    except Exception as e:
        rec["status"] = "failed"
        rec["error"] = f"{type(e).__name__}: {str(e)[:400]}"
        return rec, None


def beam_family():
    psfs = {}
    for width in BEAM_FWHM_ARCMIN:
        sigma_pix = width / (60*PIX_DEG) / (2*np.sqrt(2*np.log(2)))
        psfs[width] = gaussian_peak_psf(sigma_pix, radius_sigma=6.)
    return psfs, {
        "kind": "documented_gaussian_fallback_screen_not_measured_pa4_depth1_beam",
        "nominal_fwhm_arcmin": 1.4,
        "fwhm_sensitivity_envelope_arcmin": [1.2, 1.6],
        "envelope_is_assumed_not_measured_statistical_1sigma": True,
        "nominal_resolution_published_reference": BEAM_REFERENCE,
        "coadd_beam_archive_not_instrument_specific_depth1": BEAM_ARCHIVE_NOT_APPLICABLE,
        "psf_sha256_by_fwhm_arcmin": {str(k): digest_bytes(v.tobytes()) for k, v in psfs.items()},
        "validated_for_single_observation": False,
        "transfer_function_validated": False,
        "beam_uncertainty_propagated_into_snr_and_injection": True
    }


def time_at_position(row, ra_deg, dec_deg):
    # ACT depth-1 time FITS: ivar-weighted SECONDS AFTER PERIOD START (ctime).
    # This is a single nearest native map pixel, not a scan-duration estimate.
    path = row.get("mean_time_path")
    if not path:
        return {"status": "missing_mean_time_path"}
    try:
        if "://" in path or ".." in path or "?" in path:
            raise ValueError("Unsafe time-map path")
        url = ROOT + "/" + path.lstrip("/")
        h, meta = fits_header(url)
        x, y = WCS(h).celestial.all_world2pix(ra_deg, dec_deg, 0)
        nx, ny = int(h["NAXIS1"]), int(h["NAXIS2"])
        ix, iy = int(round(float(x))), int(round(float(y)))
        if not (0 <= ix < nx and 0 <= iy < ny):
            return {"status": "outside_time_map"}
        bpp = abs(int(h["BITPIX"]))//8
        start = int(meta["header_bytes"]) + (iy*nx+ix)*bpp
        raw, transfer = fetch_range(url, start, start+bpp-1)
        val = float(np.frombuffer(raw, dtype=">f4" if bpp == 4 else ">f8")[0])
        ctime_text = str(row["ctime"])
        if "-" in ctime_text or "T" in ctime_text or " " in ctime_text:
            period_start = datetime.fromisoformat(ctime_text).replace(tzinfo=timezone.utc).timestamp()
        else:
            period_start = float(ctime_text)
        epoch = period_start + val
        start_epoch = datetime.fromisoformat(row["start_time"]).replace(tzinfo=timezone.utc).timestamp()
        stop_epoch = datetime.fromisoformat(row["stop_time"]).replace(tzinfo=timezone.utc).timestamp()
        if not (np.isfinite(val) and 0 <= val <= stop_epoch - start_epoch + 120 and
                start_epoch-120 <= epoch <= stop_epoch+120):
            raise ValueError("Time FITS pixel disagrees with frozen scan interval")
        return {"status": "mean_time_fits_pixel_read", "unix_utc": epoch,
                "utc": datetime.fromtimestamp(epoch, timezone.utc).isoformat(),
                "period_ctime_from_catalog": row["ctime"], "period_start_unix": period_start,
                "seconds_from_period_start": val,
                "native_pixel_xy": [ix, iy], "time_map_header_sha256": meta["header_64k_sha256"],
                "time_pixel_sha256": digest_bytes(raw), "http_status": transfer["status"]}
    except Exception as e:
        return {"status": "time_pixel_failed", "error": f"{type(e).__name__}: {e}"}


def ref_coadd(ref_records):
    # Distinct-night weights: simultaneous depth-1 observations cannot be
    # treated as independent. Keep per-night image/ivar aggregates.
    nights = defaultdict(list)
    for rec, pair in ref_records:
        if pair is not None:
            nights[rec["start_time"][:10]].append(pair)
    nightly = []
    for day, maps in sorted(nights.items()):
        num = sum(m*i for m, i in maps)
        weight = sum(i for _, i in maps)
        mean = np.divide(num, weight, out=np.zeros_like(weight), where=weight>0)
        # Conservative against within-night covariance: do not sum ivars.
        ivar = np.maximum.reduce([i for _, i in maps])
        nightly.append((day, mean, ivar))
    if not nightly:
        return None, None, None
    cnt = sum((i > 0).astype(int) for _, _, i in nightly)
    denominator = sum(i for _, _, i in nightly)
    numerator = sum(m*i for _, m, i in nightly)
    out = np.divide(numerator, denominator, out=np.zeros_like(denominator), where=denominator>0)
    return out, np.where(cnt >= 3, denominator, 0), nightly


def to_float(v):
    return float(v) if v is not None and np.isfinite(v) else None



def checkpoint(result, stage):
    # If an Actions runner times out, preserve exactly which science bytes and
    # failures were observed; a partial checkpoint is not a completed search.
    result["run_stage"] = stage
    dst = BASE / "science-pilot-receipt.json"
    part = BASE / "science-pilot-receipt.json.partial"
    part.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    part.replace(dst)


def analyse(tests, refs, aperture):
    result = {"classification": "actual_fits_byte_range_pilot_not_calibrated_survey",
              "utc": datetime.now(timezone.utc).isoformat(),
              "historical_inventory_pr": 192, "new_pr_not_yet_assigned_by_runner": True,
              "eligible_gzip_sha256": MANIFEST_SHA,
              "preregistered_selection": {
                  "band": "pa4/f150", "test_start": DATE_START, "test_stop": DATE_END,
                  "reference_start": REF_START, "reference_stop": REF_END,
                  "roi_ra_deg": CENTER_RA, "roi_dec_deg": CENTER_DEC,
                  "roi_radius_deg": RADIUS_DEG, "pixel_size_arcmin": .5,
                  "target_map_ids": [str(r["map_id"]) for r in tests],
                  "reference_map_ids": [str(r["map_id"]) for r in refs]},
              "num_test_rows": len(tests), "num_reference_rows": len(refs),
              "test_observations": [], "reference_observations": [], "controls": {},
              "matched_filter": None, "injections": None,
              "effective_area_time_deg2_days": None, "rate_limit": None,
              "independent_referee_can_validate_survey": False,
              "discovery": False, "remarkable": False}
    ra, dec, region_mask = aperture
    footprint_cache = load_footprint_cache()
    result["cached_nonoverlap_header_records"] = 0
    resrefs = []
    for row in refs:
        rec, pair = read_pair(row, ra, dec, region_mask, footprint_cache)
        result["cached_nonoverlap_header_records"] += int("reuse" in rec)
        result["reference_observations"].append(rec)
        resrefs.append((rec, pair))
        print("REF", rec["map_name"], rec["status"], flush=True)
        checkpoint(result, "references_in_progress")
    reference, reference_ivar, nightly = ref_coadd(resrefs)
    result["controls"]["reference_nights_with_science_pixels"] = [day for day, _, _ in (nightly or [])]
    result["controls"]["reference_leave_out_entire_test_epoch"] = True
    testpairs = []
    for row in tests:
        rec, pair = read_pair(row, ra, dec, region_mask, footprint_cache)
        result["cached_nonoverlap_header_records"] += int("reuse" in rec)
        result["test_observations"].append(rec)
        testpairs.append((rec, pair))
        print("TEST", rec["map_name"], rec["status"], flush=True)
        checkpoint(result, "observations_in_progress")
    checkpoint(result, "before_reference_and_beam_controls")
    result["num_test_maps_science_pixels"] = sum(p is not None for _, p in testpairs)
    result["num_reference_maps_science_pixels"] = sum(p is not None for _, p in resrefs)
    if reference is None or not np.any(reference_ivar > 0):
        result["blocker"] = "No independent >=3-night reference ivar in preregistered ROI"
        return result
    psfs, beam_receipt = beam_family()
    psf = psfs[1.4]
    result["beam"] = beam_receipt
    result["controls"]["reference_epoch_leakage_count"] = sum(
        in_interval(x["start_time"], DATE_START, DATE_END) for x, p in resrefs if p is not None)
    if result["controls"]["reference_epoch_leakage_count"]:
        result["blocker"] = "Reference contamination: test epoch included"
        return result
    odd = [v for j, v in enumerate(nightly) if j % 2 == 0]
    even = [v for j, v in enumerate(nightly) if j % 2 == 1]
    split = None
    if odd and even:
        # Collapse each pre-fixed disjoint set of nights for an independent null.
        def coadd(xs):
            wi = sum(t[2] for t in xs)
            mi = np.divide(sum(t[1]*t[2] for t in xs), wi, out=np.zeros_like(wi), where=wi>0)
            return mi, wi
        m0, w0 = coadd(odd)
        m1, w1 = coadd(even)
        split = matched_filter(m0, m1, w0, w1, psf, independent_reference=True)
        result["controls"]["independent_split_night_jackknife"] = {
            "night_sets": [[r[0] for r in odd], [r[0] for r in even]],
            "negative_tail_white_7": int(np.sum(split["snr_white"] < -7)),
            "positive_tail_white_7": int(np.sum(split["snr_white"] > 7)),
            "negative_peaks_white_7": len(peaks(-np.where(split["eligible"] & region_mask, split["snr_white"], np.nan), 7, 6)),
            "positive_peaks_white_7": len(peaks(np.where(split["eligible"] & region_mask, split["snr_white"], np.nan), 7, 6)),
            "valid_centres_in_aperture": int(np.sum(split["eligible"] & region_mask)),
            "white_snr_mad": to_float(1.4826*np.nanmedian(np.abs(
                split["snr_white"][split["eligible"] & region_mask] -
                np.nanmedian(split["snr_white"][split["eligible"] & region_mask]))))
                if np.any(split["eligible"] & region_mask) else None}
    scans = []
    rng = np.random.default_rng(192)
    row_by_id = {str(row["map_id"]): row for row in tests}
    for rec, pair in testpairs:
        if pair is None:
            continue
        m, im = pair
        try:
            filtered = {k: matched_filter(m, reference, im, reference_ivar, p,
                        independent_reference=True) for k, p in psfs.items()}
        except Exception as e:
            scans.append({"map_id": rec["map_id"], "error": f"{type(e).__name__}: {e}"})
            continue
        s = filtered[1.4]
        snr = s["snr_white"]
        good = s["eligible"] & region_mask & np.isfinite(snr)
        found = peaks(np.where(good, snr, np.nan), threshold=7, distance_pix=6)
        positive = int(np.sum(good & (snr >= 7)))
        negative = int(np.sum(good & (snr <= -7)))
        candidates = []
        for p in found[:30]:
            y, x = p["y"], p["x"]
            width_screen = {}
            for width, filt in filtered.items():
                width_screen[str(width)] = {
                    "snr_white_at_nominal_peak": to_float(filt["snr_white"][y, x]),
                    "native_amplitude": to_float(filt["amplitude_native"][y, x])}
            candidates.append({**p, "ra_deg": float(ra[y,x]), "dec_deg": float(dec[y,x]),
                               "flux_native_uK_provisional": to_float(s["amplitude_native"][y,x]),
                               "gaussian_fwhm_sensitivity": width_screen,
                               "disposition": "uncalibrated_not_a_discovery_requires_catalog_and_time_vetoes"})
        yy, xx = np.mgrid[23:N-23:23, 23:N-23:23]
        choices = [(int(y), int(x)) for y, x in zip(yy.ravel(), xx.ravel()) if good[y,x]]
        if len(choices) > 25:
            choices = [choices[k] for k in sorted(rng.choice(len(choices), 25, replace=False))]
        injection_family = {}
        for width, template in psfs.items():
            filt = filtered[width]
            width_injections = {}
            for factor in (6, 8, 10, 15):
                success = 0
                for y, x in choices:
                    if not filt["eligible"][y, x]:
                        continue
                    amp = factor * float(filt["sigma_white_native"][y,x])
                    sim = inject(m, template, amp, y, x)
                    res = matched_filter(sim, reference, im, reference_ivar, template,
                                         independent_reference=True)
                    success += bool(np.isfinite(res["snr_white"][y,x]) and res["snr_white"][y,x] >= 7)
                width_injections[str(factor)] = {"recovered": success,
                                                  "attempted": len(choices)}
            injection_family[str(width)] = width_injections
        # Known Naess et al. 2021 transient: targeted, never a new discovery.
        cy = cx = N//2
        grid_y, grid_x = np.ogrid[:N, :N]
        near = ((grid_y-cy)**2 + (grid_x-cx)**2 <= 6**2) & good
        near_snr = np.where(near, snr, -np.inf)
        recovery = {"published_source": "Naess2021_2018-09-11",
                    "published_ra_deg": KNOWN_RA, "published_dec_deg": KNOWN_DEC,
                    "published_peak_time_utc_from_Naess2021": "2018-09-11T23:36:00+00:00",
                    "known_event_not_novel": True,
                    "tested_epoch_start": rec["start_time"],
                    "tested_epoch_stop": rec["stop_time"],
                    "timing": time_at_position(row_by_id[rec["map_id"]], KNOWN_RA, KNOWN_DEC)}
        if np.any(near):
            y0, x0 = np.unravel_index(int(np.argmax(near_snr)), near_snr.shape)
            separation = float(np.hypot(y0-cy,x0-cx)*0.5)
            pixel_time = recovery["timing"]
            observed_near_published_peak = (pixel_time.get("status") == "mean_time_fits_pixel_read" and
                abs(pixel_time["unix_utc"] - datetime.fromisoformat(
                    "2018-09-11T23:36:00+00:00").timestamp()) <= 6*3600)
            nominal_mad = to_float(1.4826*np.nanmedian(np.abs(
                snr[good]-np.nanmedian(snr[good])))) if np.any(good) else None
            recovery.update({"highest_white_snr_within_3arcmin": to_float(snr[y0,x0]),
                             "peak_offset_arcmin": separation,
                             "empirical_mad_rescaled_peak": to_float(
                                 snr[y0,x0]/max(nominal_mad or 1, 1)),
                             "within_2arcmin": separation <= 2.0,
                             "mean_time_within_6h_of_published_peak": bool(observed_near_published_peak),
                             "all_gaussian_fwhm_peak_snr": {
                                 str(w): to_float(f["snr_white"][y0,x0])
                                 for w,f in filtered.items()},
                             "screen_recovered_nominal": bool(snr[y0,x0] >= 7 and
                                   separation <= 2.0 and observed_near_published_peak),
                             "status": "not_a_scientific_recovery_without_independent_psf_and_noise_calibration"})
        else:
            recovery["status"] = "no_valid_filter_centres_at_published_source"
        width_tails = {
            str(w): {"positive_local_peaks_7": len(peaks(np.where(
                        f["eligible"] & region_mask, f["snr_white"], np.nan),7,6)),
                     "negative_sign_flipped_local_peaks_7": len(peaks(-np.where(
                        f["eligible"] & region_mask, f["snr_white"], np.nan),7,6))}
            for w, f in filtered.items()}
        scans.append({"map_id": rec["map_id"], "scan_time_start": rec["start_time"],
                      "scan_time_stop": rec["stop_time"],
                      "valid_psf_centres_in_aperture": int(good.sum()),
                      "positive_white_snr_pixels_ge_7": positive,
                      "sign_flip_negative_white_snr_pixels_le_minus_7": negative,
                      "positive_local_peaks_7": len(found),
                      "sign_flipped_negative_local_peaks_7": len(peaks(np.where(good, -snr, np.nan),7,6)),
                      "gaussian_beam_width_tail_sensitivity": width_tails,
                      "provisional_peaks": candidates,
                      "provisional_peaks_omitted_beyond_first_30": max(len(found)-30,0),
                      "actual_data_injections_fwhm_envelope": injection_family,
                      "published_positive_control": recovery,
                      "snr_white_mad": to_float(np.nanmedian(np.abs(snr[good] -
                          np.nanmedian(snr[good])))*1.4826) if np.any(good) else None})
    result["matched_filter"] = {
        "classification": "literature_based_gaussian_fwhm_envelope_white_noise_screen_only",
        "observation_excluded_reference": True,
        "source_noise_correlations_validated": False,
        "survey_trials_calibrated": False,
        "catalog_vetoes_run": False,
        "jpl_horizons_run": False,
        "published_transient_recovery_screen_only": any(
            s.get("published_positive_control",{}).get("screen_recovered_nominal",False) for s in scans),
        "empirical_white_snr_jackknife_mad": result["controls"].get(
            "independent_split_night_jackknife",{}).get("white_snr_mad"),
        "scans": scans}
    result["injections"] = {"kind": "actual_observation_pixels_literature_gaussian_beam_envelope_provisional",
                            "factors_white_noise_sigma": [6,8,10,15],
                            "seed": 192, "not_validated_survey_selection_function": True,
                            "fwhm_arcmin_bracket": [1.2,1.4,1.6]}
    if not scans or any("error" in s for s in scans):
        result["blocker"] = "No complete matched-filter result for at least one science-bearing test map"
    result["observed_white_noise_screened_peaks"] = sum(s.get("positive_local_peaks_7",0) for s in scans)
    result["candidate_vetting_required"] = bool(result["observed_white_noise_screened_peaks"])
    checkpoint(result, "finished_provisional_not_survey_validated")
    return result


def main():
    tests, refs = selection()
    print(f"Frozen preregistration selected {len(tests)} test + {len(refs)} reference catalog rows", flush=True)
    ra, dec, mask = tile_wcs()
    result = analyse(tests, refs, (ra, dec, mask))
    dest = BASE / "science-pilot-receipt.json"
    dest.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ("test_observations", "reference_observations", "matched_filter")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
