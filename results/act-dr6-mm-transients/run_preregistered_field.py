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
from mf_core import matched_filter, peaks, inject

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
# Published DR6 coadd beam (measured), not validated for a single depth-1
# observation; any corresponding MF SNR is PROVISIONAL only.
BEAM_URL = "https://lambda.gsfc.nasa.gov/data/act/beams/act_dr6.02_daytime_beams.tar.gz"


def selection():
    if digest(MANIFEST) != MANIFEST_SHA:
        raise ValueError("Frozen compressed catalog hash mismatch")
    rows = [json.loads(line) for line in gzip.open(MANIFEST, "rt", encoding="utf-8")]
    rows = [r for r in rows if r["tube_slot"] == "pa4" and r["frequency"] == "f150"]
    tests = [r for r in rows if DATE_START <= r["start_time"] < DATE_END]
    refs = [r for r in rows if REF_START <= r["start_time"] < REF_END
            and not (DATE_START <= r["start_time"] < DATE_END)]
    if len({r["map_id"] for r in tests + refs}) != len(tests + refs):
        raise ValueError("Duplicate map id in frozen subset")
    return tests, refs


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


def access_tile(path, ra, dec):
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


def read_pair(row, ra, dec, aperture):
    rec = {"map_id": str(row["map_id"]), "map_name": row["map_name"],
           "start_time": row["start_time"], "stop_time": row["stop_time"]}
    try:
        m, rm = access_tile(row["map_path"], ra, dec)
        i, ri = access_tile(row["ivar_path"], ra, dec)
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


def measured_coadd_psf():
    """Attempt public measured beam profile; not depth-1 beam validation."""
    req = Request(BEAM_URL, headers={"User-Agent": USER_AGENT})
    with urlopen(req, timeout=90) as response:
        data = response.read(25*1024*1024 + 1)
    if len(data) > 25*1024*1024:
        raise ValueError("Beam archive exceeded 25 MiB")
    names = []
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tar:
        for item in tar:
            if item.isfile() and "pa4_f150" in item.name and "coadd" in item.name and ("jitter_mono" in item.name or "instant" in item.name):
                names.append(item.name)
        if not names:
            raise ValueError("pa4_f150 measured coadd profile missing from beam archive")
        name = sorted(names, key=lambda n: ("jitter_mono" not in n, len(n), n))[0]
        raw = tar.extractfile(name).read()
    arr = np.loadtxt(io.BytesIO(raw), comments="#")
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError("Beam profile not two-column")
    ell, bl = arr[:, 0], arr[:, 1]
    order = np.argsort(ell)
    ell, bl = ell[order], bl[order]
    keep = (ell >= 0) & (ell <= 30000) & np.isfinite(bl)
    ell, bl = ell[keep], bl[keep]
    if len(ell) < 50:
        raise ValueError("Insufficient measured beam multipoles")
    yy, xx = np.mgrid[-12:13, -12:13]
    theta = np.hypot(xx, yy)*np.deg2rad(PIX_DEG)
    profile = np.trapezoid(ell[:,None,None]*bl[:,None,None]*j0(ell[:,None,None]*theta[None,:,:]), ell, axis=0)
    profile /= profile[12, 12]
    min_value = float(profile.min())
    # The white-noise pilot uses the nonnegative central lobe only.
    psf = np.maximum(profile, 0.)
    psf /= psf.max()
    return psf, {"kind": "published_DR6_coadd_harmonic_beam_not_validated_for_depth1",
                 "beam_archive_url": BEAM_URL, "archive_sha256": digest_bytes(data),
                 "beam_member": name, "member_sha256": digest_bytes(raw),
                 "negative_sidelobe_min_before_clipping": min_value,
                 "psf_sha256": digest_bytes(psf.tobytes()), "validated_for_single_observation": False}


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
    part.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\\n")
    part.replace(dst)


def analyse(tests, refs, aperture):
    result = {"classification": "actual_fits_byte_range_pilot_not_calibrated_survey",
              "utc": datetime.now(timezone.utc).isoformat(),
              "pr_number": 192, "pr_already_merged_before_this_experiment": True,
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
    resrefs = []
    for row in refs:
        rec, pair = read_pair(row, ra, dec, region_mask)
        result["reference_observations"].append(rec)
        resrefs.append((rec, pair))
        print("REF", rec["map_name"], rec["status"], flush=True)
        checkpoint(result, "references_in_progress")
    reference, reference_ivar, nightly = ref_coadd(resrefs)
    result["controls"]["reference_nights_with_science_pixels"] = [day for day, _, _ in (nightly or [])]
    result["controls"]["reference_leave_out_entire_test_epoch"] = True
    testpairs = []
    for row in tests:
        rec, pair = read_pair(row, ra, dec, region_mask)
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
    try:
        psf, beam_receipt = measured_coadd_psf()
        result["beam"] = beam_receipt
    except Exception as e:
        result["beam"] = {"status": "measured_beam_unavailable", "error": f"{type(e).__name__}: {e}"}
        result["blocker"] = "Published measured beam profile not retrieved; no filter run"
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
            "valid_centres": int(split["eligible"].sum())}
    scans = []
    rng = np.random.default_rng(192)
    for rec, pair in testpairs:
        if pair is None:
            continue
        m, im = pair
        try:
            s = matched_filter(m, reference, im, reference_ivar, psf, independent_reference=True)
        except Exception as e:
            scans.append({"map_id": rec["map_id"], "error": f"{type(e).__name__}: {e}"})
            continue
        snr = s["snr_white"]
        good = s["eligible"] & region_mask & np.isfinite(snr)
        found = peaks(np.where(good, snr, np.nan), threshold=7, distance_pix=6)
        positive = int(np.sum(good & (snr >= 7)))
        negative = int(np.sum(good & (snr <= -7)))
        candidates = []
        for p in found[:30]:
            y, x = p["y"], p["x"]
            candidates.append({**p, "ra_deg": float(ra[y,x]), "dec_deg": float(dec[y,x]),
                               "flux_native_uK_provisional": to_float(s["amplitude_native"][y,x]),
                               "disposition": "unvetted_noise_uncalibrated_peak_not_discovery"})
        yy, xx = np.mgrid[23:N-23:23, 23:N-23:23]
        choices = [(int(y), int(x)) for y, x in zip(yy.ravel(), xx.ravel()) if good[y,x]]
        if len(choices) > 25:
            choices = [choices[k] for k in sorted(rng.choice(len(choices), 25, replace=False))]
        inj = {}
        for factor in (6, 8, 10, 15):
            success = 0
            for y, x in choices:
                amp = factor * float(s["sigma_white_native"][y,x])
                test = inject(m, psf, amp, y, x)
                res = matched_filter(test, reference, im, reference_ivar, psf, independent_reference=True)
                success += bool(np.isfinite(res["snr_white"][y,x]) and res["snr_white"][y,x] >= 7)
            inj[str(factor)] = {"recovered": int(success), "attempted": len(choices)}
        scans.append({"map_id": rec["map_id"], "valid_psf_centres_in_aperture": int(good.sum()),
                      "positive_white_snr_pixels_ge_7": positive,
                      "sign_flip_negative_white_snr_pixels_le_minus_7": negative,
                      "provisional_peaks": candidates, "actual_data_injections": inj,
                      "snr_white_mad": to_float(np.nanmedian(np.abs(snr[good] -
                          np.nanmedian(snr[good])))*1.4826) if np.any(good) else None})
    result["matched_filter"] = {
        "classification": "published_DR6_coadd_beam_white_noise_screen_only",
        "observation_excluded_reference": True,
        "source_noise_correlations_validated": False,
        "survey_trials_calibrated": False,
        "catalog_vetoes_run": False,
        "jpl_horizons_run": False,
        "published_transient_recovered": False,
        "scans": scans}
    result["injections"] = {"kind": "actual_observation_pixels_coadd_beam_provisional",
                            "factors_white_noise_sigma": [6,8,10,15],
                            "seed": 192, "not_validated_survey_selection_function": True}
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
