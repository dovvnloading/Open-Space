#!/usr/bin/env python3
"""Control-first XMM QPE validation. No archive search unless controls pass."""
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import tarfile
import tempfile
import time
from urllib.parse import urlencode

import numpy as np
import requests
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)
MANIFEST = OUT / "frozen_selection.json"
# Only XMM epochs explicitly reported as QPE-active by the cited publications.
# AT2019qiz's confirmed QPEs were observed by NICER/Swift/AstroSat, not XMM.
# Ansky's >4-day eruptions were reported after the DR14 release.
TARGETS = {
    "GSN 069": {
        "obsid": "0851180401", "published_h": [7.0, 11.0],
        "source": "https://doi.org/10.1051/0004-6361/202244512",
        "note": "2019 May 31, QPE-active (Kosec et al. Table 4)"},
    "RX J1301.9+2747": {
        "obsid": "0864560101", "published_h": [1.0, 6.0],
        "source": "https://doi.org/10.1051/0004-6361/202450861",
        "note": "Irregular QPE peaks, ~4 h mean"},
    "eRO-QPE1": {
        "obsid": "0861910201", "published_h": [1.5, 30.0],
        "source": "https://doi.org/10.1038/s41586-021-03394-6",
        "note": "2020 July 27 XMM1, three overlapping eruptions; do not equate subpeaks with period"},
    "eRO-QPE2": {
        "obsid": "0872390101", "published_h": [1.6, 4.0],
        "source": "https://doi.org/10.1051/0004-6361/202349234",
        "note": "2020 Aug QPEs; typical ~2.4 h, not a strict clock"},
    "eRO-QPE3": {
        "obsid": "0883770101", "published_h": [16.0, 22.0],
        "source": "https://doi.org/10.1051/0004-6361/202348881",
        "note": "2022 July XMM1, two peaks ~20.4 h apart"},
    "eRO-QPE4": {
        "obsid": "0883770401", "published_h": [9.0, 16.0],
        "source": "https://doi.org/10.1051/0004-6361/202348881",
        "note": "2023 March XMM, repeat ~9.8-15.5 h; PPS identity needs checking"},
}
# Published X-ray astrometry, not pointing coordinates. Header SRC_RA/SRC_DEC is independently checked.
POS = {
    "GSN 069": (19.786098, -34.191813),
    "RX J1301.9+2747": (195.50054, 27.78267),
    "eRO-QPE1": (37.946917, -10.336197),
    "eRO-QPE2": (38.704042, -44.325458),
    "eRO-QPE3": (210.2222, -28.7665),
    "eRO-QPE4": (71.3921, -10.2005),  # Arcodia et al. 2024, X-ray position
}

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "Open-Space-XMM-QPE-validated-controls/2.0"})
MAX_TAR_BYTES = 3 * 1024**3  # Not the failed 120 MiB cap.
MAX_SECONDS_PER_TARGET = 14 * 60  # Three targets plus installation remain within 55min
BINSIZE_S = 300
RECEIPTS = []


def safe(obj):
    """Strict JSON: nonfinite FITS catalogue entries become JSON null."""
    if isinstance(obj, dict):
        return {k: safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [safe(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return float(obj) if math.isfinite(float(obj)) else None
    if isinstance(obj, np.integer):
        return int(obj)
    return obj


def store(name, data):
    (OUT / name).write_text(json.dumps(safe(data), indent=2, sort_keys=True, allow_nan=False) + "\n")


def checksum(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024**2), b""):
            h.update(block)
    return h.hexdigest()


def download_product(obsid, target, workdir):
    params = {"obsno": obsid, "level": "PPS", "extension": "FTZ",
              "instname": "PN", "name": "SRCTSR"}
    url = "https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?" + urlencode(params)
    dest = workdir / "source_series.tar"
    rec = {"control": target, "obsid": obsid, "url": url,
           "requested_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "byte_cap": MAX_TAR_BYTES, "wall_cap_s": MAX_SECONDS_PER_TARGET}
    t0 = time.monotonic()
    h = hashlib.sha256()
    try:
        with SESSION.get(url, stream=True, timeout=(25, 100)) as response:
            rec.update(http_status=response.status_code,
                       final_url=response.url,
                       content_length=response.headers.get("Content-Length"),
                       content_type=response.headers.get("Content-Type"))
            response.raise_for_status()
            n = 0
            if int(response.headers.get("Content-Length") or 0) > MAX_TAR_BYTES:
                raise ValueError("PPS archive over 3 GiB safety ceiling")
            with dest.open("wb") as f:
                for chunk in response.iter_content(1024**2):
                    if not chunk:
                        continue
                    n += len(chunk)
                    if n > MAX_TAR_BYTES:
                        raise ValueError("PPS stream over 3 GiB safety ceiling")
                    if time.monotonic() - t0 > MAX_SECONDS_PER_TARGET:
                        raise TimeoutError("control request exceeded wall-clock allocation")
                    h.update(chunk)
                    f.write(chunk)
            if n < 1024:
                raise ValueError("implausibly small tar file")
            rec.update(bytes=n, sha256=h.hexdigest(), status="downloaded")
            return dest
    except Exception as exc:
        rec.update(status="unavailable", error=repr(exc)[:400])
        dest.unlink(missing_ok=True)
        return None
    finally:
        rec["elapsed_s"] = round(time.monotonic() - t0, 1)
        RECEIPTS.append(rec)


def get_manifest():
    # json.loads accepts the legacy NaN tokens; store() repairs them to null.
    with MANIFEST.open() as f:
        m = json.load(f)
    store("frozen_selection.json", m)
    return safe(m)


def identify(target, manifest):
    """Use original DR14 position-matched source number; no guessed PPS identity."""
    spec = TARGETS[target]
    rows = manifest.get("known_qpe_matches", {}).get(target, [])
    for row in rows:
        if str(row["obs_id"]).zfill(10) == spec["obsid"]:
            return {"pps_srcnum": int(row["pps_srcnum"]),
                    "catalog_ra": row.get("ra"), "catalog_dec": row.get("dec"),
                    "detid": row.get("detid"),
                    "identity": "DR14 position-matched within 12 arcsec; PPS association needs header verification"}
    return {"pps_srcnum": None, "identity": "NO DR14 positional match for specified QPE-active observation"}


def mjd_reference(header):
    if "MJDREF" in header:
        return float(header["MJDREF"])
    return float(header.get("MJDREFI", 50814.0)) + float(header.get("MJDREFF", 0))


def combine_energies(hdu):
    """PNX000 source RATE1-3 are net 0.2-2 keV, BACK1-3 already scaled."""
    d, keys = hdu.data, set(hdu.data.names or [])
    needed = ["TIME", "FRACEXP"] + [f"RATE{i}" for i in (1, 2, 3)] + [
        f"RATE{i}_ERR" for i in (1, 2, 3)] + [
        f"BACK{i}V" for i in (1, 2, 3)]
    missing = [c for c in needed if c not in keys]
    if missing:
        raise ValueError("No band-resolved merged PN net series; missing " + ",".join(missing))
    net = sum(np.asarray(d[f"RATE{i}"], dtype=float) for i in (1, 2, 3))
    variance = sum(np.asarray(d[f"RATE{i}_ERR"], dtype=float)**2 for i in (1, 2, 3))
    bkg = sum(np.asarray(d[f"BACK{i}V"], dtype=float) for i in (1, 2, 3))
    hard = None
    if "RATE4" in keys and "RATE5" in keys:
        hard = np.asarray(d["RATE4"], float) + np.asarray(d["RATE5"], float)
    timedel = np.asarray(d["TIMEDEL"], float) if "TIMEDEL" in keys else None
    return np.asarray(d["TIME"], float), net, np.sqrt(np.maximum(variance, 0)), bkg, (
        np.asarray(d["FRACEXP"], float)), hard, timedel


def rebin(time_s, rate, error, bkg, fracexp, hard=None, timedel=None, step=BINSIZE_S):
    """Exposure-weighted rate/error aggregation; no interpolation across GTI gaps."""
    mask = (np.isfinite(time_s) & np.isfinite(rate) & np.isfinite(error) &
            (error >= 0) & np.isfinite(bkg) & np.isfinite(fracexp) & (fracexp >= .4))
    if np.count_nonzero(mask) < 100:
        raise ValueError("Insufficient valid native EPIC bins")
    t, y, e, b, f = [a[mask] for a in (time_s, rate, error, bkg, fracexp)]
    order = np.argsort(t)
    t, y, e, b, f = [a[order] for a in (t, y, e, b, f)]
    dtarr = np.diff(t)
    native_dt = np.nanmedian(dtarr[(dtarr > 0) & (dtarr < 100)])
    if not np.isfinite(native_dt) or native_dt <= 0:
        raise ValueError("Invalid native bin step")
    if timedel is not None:
        td = np.asarray(timedel[mask], float)[order]
        td = np.where(np.isfinite(td) & (td > 0) & (td < 1000), td, native_dt)
    else:
        td = np.full(len(t), native_dt)
    exp = td * f
    start = math.floor(float(t[0]) / step) * step
    bins = start + np.arange(int(math.ceil((t[-1] - start) / step)) + 2) * step
    norm = np.histogram(t, bins=bins, weights=exp)[0]
    val = np.histogram(t, bins=bins, weights=y * exp)[0]
    err2 = np.histogram(t, bins=bins, weights=(e * exp)**2)[0]
    bgval = np.histogram(t, bins=bins, weights=b * exp)[0]
    good = norm >= .65 * step
    with np.errstate(divide="ignore", invalid="ignore"):
        rr = np.where(good, val / norm, np.nan)
        ee = np.where(good, np.sqrt(err2) / norm, np.nan)
        bb = np.where(good, bgval / norm, np.nan)
    hh = None
    if hard is not None:
        hv = np.asarray(hard[mask], float)[order]
        hv = np.where(np.isfinite(hv), hv, 0)
        hval = np.histogram(t, bins=bins, weights=hv * exp)[0]
        hh = np.where(good, hval / np.maximum(norm, 1), np.nan)
    return {"t_s": bins[:-1] + .5 * step, "rate": rr, "error": ee,
            "background": bb, "hard": hh, "exposure": norm, "good": good,
            "bin_s": step, "native_bin_s": round(float(native_dt), 3)}


def find_eruptions(curve, refmjd, allowed_period):
    """Control-calibrated 300s matched-filter screen, NOT independent QPE confirmation."""
    t, y, e, bg, good = [curve[k] for k in ("t_s", "rate", "error", "background", "good")]
    idx = np.flatnonzero(good & np.isfinite(y) & np.isfinite(e) & np.isfinite(bg))
    if len(idx) < 24:
        return {"status": "insufficient_exposure", "good_bins": len(idx)}
    yy = y[idx]
    quiet = float(np.percentile(yy, 20))
    bquiet = float(np.percentile(bg[idx], 30))
    errtyp = max(float(np.median(e[idx])), 1e-5)
    mad = float(1.4826 * np.median(abs(yy - np.median(yy))))
    # Smooth only within good-time exposure; downweight missing data and reject edge candidates.
    w = good.astype(float)
    smooth = gaussian_filter1d(np.nan_to_num(y) * w, 1.2, mode="constant") / np.maximum(
        gaussian_filter1d(w, 1.2, mode="constant"), 1e-8)
    smooth[w < .5] = np.nan
    peak_min = quiet + max(4 * errtyp, 2.0 * mad)
    prominence = max(4 * errtyp, 1.0 * mad)
    candidate_idx, properties = find_peaks(
        np.nan_to_num(smooth, nan=-1e9), height=peak_min,
        prominence=prominence, distance=12)  # >=60 min
    accepted = []
    for p in candidate_idx:
        lo, hi = max(0, p-2), min(len(t), p+3)
        if p < 2 or p + 2 >= len(t) or not all(good[lo:hi]):
            continue
        local = slice(max(0, p-4), min(len(t), p+5))
        net_excess = float(smooth[p] - quiet)
        bgpeak = float(np.nanmedian(bg[local]))
        bge = bgpeak - bquiet
        snr = net_excess / max(float(e[p]), errtyp/2, 1e-6)
        # Broad flares span many 300s bins; using single-bin >=5 sigma
        # missed the second published eRO-QPE3 eruption despite a clear plateau.
        wix = np.arange(max(0, p-6), min(len(t), p+7))
        wix = wix[good[wix] & np.isfinite(y[wix]) & np.isfinite(e[wix])]
        kernel = np.exp(-.5*((t[wix]-t[p])/850.)**2)
        variance = np.maximum(e[wix], max(errtyp/2, 1e-6))**2
        integrated = float(np.sum(kernel*(y[wix]-quiet)/variance) /
                           np.sqrt(np.sum(kernel**2/variance))) if len(wix)>=6 else 0.
        ratio = max(float(smooth[p]), 0) / max(quiet, .001)
        background_veto = bge > max(.5*net_excess, 5 * np.nanmedian(e[idx]))
        if snr < 3.5 or integrated < 7 or ratio < 3 or background_veto:
            continue
        accepted.append({"mjd_tt": round(float(refmjd + t[p]/86400), 6),
                         "time_since_first_h": round(float((t[p]-t[idx[0]])/3600), 3),
                         "net_rate_cps": round(float(smooth[p]), 5),
                         "excess_sigma": round(snr, 2),
                         "integrated_sigma": round(integrated, 2),
                         "amplitude_vs_quiet": round(ratio, 2),
                         "background_rate_cps": round(bgpeak, 5)})
    dts = np.diff([p["mjd_tt"] for p in accepted]) * 24
    interval = [round(float(x), 3) for x in dts]
    match = len(dts) > 0 and any(allowed_period[0] <= x <= allowed_period[1] for x in dts)
    return {"status": "recurrence_range_match_preliminary" if match else "no_matched_recurring_train",
            "n_peaks": len(accepted), "peaks": accepted,
            "intervals_h": interval, "published_range_h": allowed_period,
            "range_match": bool(match), "net_baseline_cps": round(quiet, 6),
            "background_baseline_cps": round(bquiet, 6),
            "binned_good_exposure_h": round(float(np.sum(curve["exposure"][good])/3600), 3),
            "span_h": round(float((t[idx[-1]]-t[idx[0]])/3600), 3),
            "native_bin_s": curve["native_bin_s"], "bin_s": curve["bin_s"],
            "background_veto": "source BACK1-3 compared locally; FBKTSR not yet independently checked",
            "published_peak_time_match": "NOT VERIFIED (ranges are not individual timestamps)",
            "pile_up_status": "NOT VALIDATED",
            "spectral_temperature_evolution": "NOT MEASURED (energy-resolved count rates only)"}


def save_binned_curve(target, obsid, c, refmjd):
    path = OUT / f"control_curve_{target.replace(' ', '_').replace('+', 'p')}_{obsid}.csv"
    with path.open("w", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(["mjd_tt", "time_s", "net_0p2_2keV_cps", "net_err_cps",
                         "scaled_background_0p2_2keV_cps", "hard_2_12keV_cps", "exposure_s"])
        for i, t in enumerate(c["t_s"]):
            if not c["good"][i]:
                continue
            writer.writerow([round(float(refmjd + t/86400), 8),
                             round(float(t), 1), round(float(c["rate"][i]), 7),
                             round(float(c["error"][i]), 7),
                             round(float(c["background"][i]), 7),
                             round(float(c["hard"][i]), 7) if c["hard"] is not None else "",
                             round(float(c["exposure"][i]), 2)])
    return {"path": str(path.relative_to(BASE)), "sha256": checksum(path)}


def source_sky_position(headers):
    """The observation pointing is not the position of a detected source."""
    for header in headers:
        for ra_key, dec_key in (("SRC_RA", "SRC_DEC"), ("RA_OBJ", "DEC_OBJ"),
                                ("RA", "DEC")):
            if ra_key not in header or dec_key not in header:
                continue
            try:
                ra, dec = float(header[ra_key]), float(header[dec_key])
                if 0 <= ra < 360 and -90 <= dec <= 90:
                    return ra, dec, ra_key, dec_key
            except (TypeError, ValueError):
                pass
    return None


def sep_arcsec(pos, expected):
    ra1, dec1 = pos[:2]
    ra2, dec2 = expected
    delta_ra = ((ra1 - ra2 + 180) % 360 - 180) * math.cos(math.radians((dec1+dec2)/2))
    return 3600 * math.hypot(delta_ra, dec1-dec2)


def analyze_tar(tarpath, target, manifest):
    spec = TARGETS[target]
    ident = identify(target, manifest)
    result = {"control": target, "obsid": spec["obsid"], **ident,
              "published_period_range_h": spec["published_h"],
              "published_reference": spec["source"], "published_note": spec["note"],
              "expected_ra_deg": POS[target][0], "expected_dec_deg": POS[target][1],
              "status": "not_analyzed"}
    found = []
    # Select the merged observation-wide band-resolved PPS source series by actual
    # source coordinates, because the latest XSA PPS source indexes may differ from DR14.
    # A headerless file may be inspected diagnostically, never taken as a verified match.
    with tarfile.open(tarpath, "r:*") as archive:
        candidates = [m for m in archive if m.isfile() and m.size <= 100*1024**2 and
                      Path(m.name).name.upper().startswith("P" + spec["obsid"] + "PNX000SRCTSR8")
                      and Path(m.name).suffix.upper() in (".FTZ", ".FIT", ".FITS")]
        result["n_merged_pn_products"] = len(candidates)
        result["available_merged_products"] = [
            {"name": Path(m.name).name, "compressed_size": m.size} for m in candidates[:50]]
        for m in candidates:
            name = Path(m.name).name.upper()
            expected_idx = f"{ident['pps_srcnum']:03X}" if ident["pps_srcnum"] is not None else None
            from_catalog = expected_idx is not None and expected_idx in name[-7:]
            try:
                zipped = archive.extractfile(m).read()
                raw = gzip.decompress(zipped) if name.endswith(".FTZ") else zipped
                with fits.open(io.BytesIO(raw), memmap=False) as hdul:
                    rates = next((x for x in hdul[1:] if x.data is not None and
                                  hasattr(x.data, "names") and "RATE1" in x.data.names), None)
                    hdrs = [rates.header, hdul[0].header] if rates is not None else [hdul[0].header]
                    pos = source_sky_position(hdrs)
                    separation = sep_arcsec(pos, POS[target]) if pos else None
                    info = {"name": name, "index_matches_DR14": bool(from_catalog),
                            "source_position": list(pos[:2]) if pos else None,
                            "position_keywords": list(pos[2:]) if pos else None,
                            "separation_arcsec": round(separation, 2) if separation is not None else None,
                            "has_energy_bands": bool(rates)}
                    found.append(info)
                    # Exact matched source position is decisive. Do not silently trust
                    # a DR14 source index if the current PPS source position differs.
                    if separation is None or separation > 15.0 or rates is None:
                        continue
                    if str(rates.header.get("HDUCLAS2", "")).strip() != "NET":
                        raise ValueError("PPS source rate not certified background-subtracted NET")
                    t, net, err, bg, frac, hard, timedel = combine_energies(rates)
                    refmjd = mjd_reference(rates.header)
                    curve = rebin(t, net, err, bg, frac, hard, timedel)
                    result["analysis"] = find_eruptions(curve, refmjd, spec["published_h"])
                    result["lightcurve"] = save_binned_curve(target, spec["obsid"], curve, refmjd)
                    result.update(product_sha256=hashlib.sha256(raw).hexdigest(),
                                  product_filename=name, columns=list(rates.data.names),
                                  hduclass2=rates.header.get("HDUCLAS2"), mjdref_tt=refmjd,
                                  status="processed_position_verified_net_pn",
                                  source_position_deg=list(pos[:2]),
                                  source_separation_arcsec=round(separation, 2),
                                  identity="PPS RATE/primary FITS source coordinates matched published X-ray position")
                    break
            except Exception as exc:
                found.append({"name": name, "parse_error": repr(exc)[:250]})
    result["product_header_audit"] = found[:60]
    if result["status"] == "not_analyzed":
        result["status"] = (
            "no_merged_pn_series_in_pps" if not candidates else
            "no_fits_position_verified_soft_source_series")
    return result

def verify_synthetic():
    """Do not run archive if detector fails bright, dim or background-only fixtures."""
    tt = np.arange(0, 38*3600, BINSIZE_S, dtype=float)
    base = np.full(tt.size, .035)
    bb = np.full(tt.size, .015)
    e = np.full(tt.size, .012)
    flares = sum(.9 * np.exp(-.5*((tt-(h*3600))/(1200))**2) for h in (4, 14, 24, 34))
    c = {"t_s": tt, "rate": base+flares, "error": e, "background": bb,
         "good": np.ones(tt.size, dtype=bool), "exposure": np.full(tt.size, BINSIZE_S),
         "bin_s": BINSIZE_S, "native_bin_s": 2.0}
    a = find_eruptions(c, 58000, (8, 12))
    if a.get("n_peaks") != 4 or not a.get("range_match"):
        raise AssertionError("synthetic source QPE control failed: " + repr(a))
    c["background"] = bb + flares
    b = find_eruptions(c, 58000, (8, 12))
    if b.get("n_peaks", 0) != 0:
        raise AssertionError("background flare veto failed: " + repr(b))
    return {"four_peak_10h_fixture": True, "background_only_fixture_rejected": True}


def main():
    batch = json.loads((BASE / "control_batch.json").read_text())
    if batch.get("stage") != "controls":
        raise SystemExit("Archive gate is closed; set stage=controls for control recovery first.")
    targets = batch["targets"]
    assert len(targets) <= 3 and len(set(targets)) == len(targets), "Three controls max per 55min job"
    tag = batch["batch_id"]
    assert tag.isidentifier() and all(t in TARGETS for t in targets)
    test = verify_synthetic()
    manifest = get_manifest()
    rows = []
    for target in targets:
        spec = TARGETS[target]
        print("CONTROL", target, spec["obsid"], flush=True)
        with tempfile.TemporaryDirectory() as td:
            tarpath = download_product(spec["obsid"], target, Path(td))
            if tarpath is None:
                rows.append({"control": target, "obsid": spec["obsid"], "status": "PPS_download_failed",
                             "error": RECEIPTS[-1].get("error")})
            else:
                try:
                    rows.append(analyze_tar(tarpath, target, manifest))
                except Exception as exc:
                    rows.append({"control": target, "obsid": spec["obsid"],
                                 "status": "tar_analysis_error", "error": repr(exc)})
        store(f"controls_{tag}.json", {"method": "NET PN merged RATE1+2+3, 0.2-2 keV, 300s",
                                     "synthetic": test, "controls": rows})
        store(f"control_receipts_{tag}.json", RECEIPTS)
        print("CONTROL RESULT", target, rows[-1].get("status"),
              rows[-1].get("analysis", {}).get("n_peaks"), flush=True)
    # Do not claim per-control confirmed recovery from broad recurrence ranges alone.
    store(f"control_gate_{tag}.json",
          {"controls_processed": len(rows),
           "recurrence_range_matches": sum(
               bool(r.get("analysis", {}).get("range_match")) for r in rows),
           "published_peak_times_independently_verified": 0,
           "archive_screen_permitted": False,
           "strict_json": True, "schema": "control-first-v2"})


if __name__ == "__main__":
    main()
