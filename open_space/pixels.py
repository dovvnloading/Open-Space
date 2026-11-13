"""Localize a periodic dimming signal in calibrated TESS SPOC postage stamps."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
from urllib.request import urlopen
import warnings

import numpy as np


PRODUCT = re.compile(r"tess\d{13}-s(\d{4})-(\d{16})-\d{4}-s_(?:lc|tp)\.fits\Z")
MAX_BYTES = 400 * 1024 * 1024
PHASE_OFFSETS = (-0.25, 0.25)
METHOD = "tess-bracketed-difference-v1"
TRANSIT_WINDOWS = ((-.4, .4), (-1.8, -.8), (.8, 1.8))
VARIABLE_WINDOWS = ((-.10, .10), (-.50, -.40), (.40, .50))


def _digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate(task):
    name = task["product_name"]
    match = PRODUCT.fullmatch(name)
    if not match or int(match[1]) != task["sector"] or int(match[2]) != task["tic_id"]:
        raise ValueError("SPOC basename, TIC ID and sector must agree")
    period, epoch, duration = (float(task[k]) for k in ("period_days", "epoch_btjd", "duration_days"))
    if not all(map(math.isfinite, (period, epoch, duration))) or period <= 0 or not 0 < duration < period / 4:
        raise ValueError("Supply a finite BTJD epoch, a positive period and a duration below one quarter period")
    expected = task.get("source_sha256")
    if expected is not None and not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ValueError("source_sha256 must be the lower-case SHA-256 for the target pixel FITS")
    return name.replace("_lc.fits", "_tp.fits"), period, epoch, duration


def _acquire(cache, name, expected):
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / name
    uri = "mast:TESS/product/" + name
    url = "https://mast.stsci.edu/api/v0.1/Download/file?uri=" + uri
    if not path.exists():
        td = tempfile.NamedTemporaryFile(prefix=".tp-", dir=cache, delete=False)
        temp = Path(td.name)
        try:
            with td, urlopen(url, timeout=90) as response:
                total = 0
                while block := response.read(1024 * 1024):
                    total += len(block)
                    if total > MAX_BYTES:
                        raise ValueError("Target pixel file exceeds the per-product memory limit")
                    td.write(block)
            if expected and expected != _digest(temp):
                raise ValueError("Downloaded target pixel checksum mismatch")
            # Concurrent workers race only to publish complete immutable products.
            try:
                os.link(temp, path)
            except FileExistsError:
                if _digest(temp) != _digest(path):
                    raise ValueError("Conflicting bytes at the same MAST product name")
        finally:
            temp.unlink(missing_ok=True)
    if path.stat().st_size > MAX_BYTES:
        raise ValueError("Cached target pixel file exceeds the per-product memory limit")
    sha = _digest(path)
    if expected and expected != sha:
        raise ValueError("Cached target pixel checksum mismatch")
    return path, uri, url, sha


def _windows(times, epoch, period, duration, cadence, intervals=TRANSIT_WINDOWS):
    """Select complete two-sided brackets at pre-specified image-independent epochs."""
    output = []
    start = math.ceil((times[0] - epoch) / period)
    stop = math.floor((times[-1] - epoch) / period)
    if stop - start > len(times):
        raise ValueError("Event frequency exceeds cadence resolution")
    for index in range(start, stop + 1):
        center = epoch + index * period
        # Use the central 80% and separated sidebands, avoiding ingress/egress.
        windows = []
        for lo, hi in intervals:
            bounds = center + duration * np.array((lo, hi))
            a, b = np.searchsorted(times, bounds)
            required = max(3, math.ceil(0.55 * (hi-lo) * duration / cadence))
            if b - a < required:
                break
            windows.append(slice(a, b))
        if len(windows) == 3:
            output.append((index, center, windows))
    return output


def _estimated_image(data, errors, sl):
    cube, uncertainty = data[sl], errors[sl]
    valid = np.isfinite(cube) & np.isfinite(uncertainty) & (uncertainty > 0)
    n = valid.sum(axis=0)
    enough = n >= max(3, math.ceil(0.75 * len(cube)))
    mean = np.full(cube.shape[1:], np.nan)
    var = mean.copy()
    mean[enough] = np.where(valid, cube, 0).sum(axis=0, dtype=float)[enough] / n[enough]
    var[enough] = np.where(valid, uncertainty**2, 0).sum(axis=0, dtype=float)[enough] / n[enough]**2
    return mean, var


def event_images(t, cube, err, epoch, period, duration, cadence, intervals=TRANSIT_WINDOWS):
    """Measure flux deficits and propagate per-cadence errors in e-/s squared."""
    records, images, variances, refs = [], [], [], []
    for cycle, center, windows in _windows(t, epoch, period, duration, cadence, intervals):
        inside, left, right = (_estimated_image(cube, err, sl) for sl in windows)
        tm, tl, tr = (np.mean(t[sl]) for sl in windows)
        wr = (tm - tl) / (tr - tl)
        out = (1 - wr)*left[0] + wr*right[0]
        diff = out - inside[0]
        var = inside[1] + (1-wr)**2 * left[1] + wr**2 * right[1]
        records.append({"cycle": cycle, "center_btjd": float(center),
                        "cadences_in_left_right": [int(s.stop - s.start) for s in windows]})
        images.append(diff)
        variances.append(var)
        refs.append(out)
    return records, images, variances, refs


def _combine(stack, variances):
    weights = np.zeros_like(variances, dtype=float)
    valid = np.isfinite(stack) & np.isfinite(variances) & (variances > 0)
    np.divide(1.0, variances, where=valid, out=weights)
    total = weights.sum(axis=0)
    output = np.full(stack.shape[1:], np.nan)
    # Require at least half of the accepted events to contribute to this pixel.
    use = (total > 0) & (valid.sum(axis=0) >= math.ceil(len(stack)/2))
    output[use] = (np.where(valid, stack, 0) * weights).sum(axis=0)[use] / total[use]
    photon_se = np.full_like(output, np.nan)
    photon_se[use] = 1 / np.sqrt(total[use])
    if len(stack) > 1:
        scatter = np.nansum((np.where(valid, stack, np.nan) - output)**2 * weights**2, axis=0)
        count = valid.sum(axis=0)
        err = np.sqrt(scatter * count / np.maximum(count - 1, 1)) / np.where(use, total, np.nan)
        photon_se = np.fmax(photon_se, err)
    return output, photon_se


def _centroid(image, peak=None):
    valid = np.isfinite(image)
    if not valid.any():
        return None
    level = float(np.nanmedian(image))
    signal = np.where(valid, image - level, -np.inf)
    if peak is None:
        peak = np.unravel_index(np.argmax(signal), signal.shape)
    y, x = np.indices(image.shape)
    selected = (abs(x - peak[1]) <= 1) & (abs(y - peak[0]) <= 1) & (signal > 0)
    weights = np.where(selected, signal, 0)
    total = weights.sum()
    if total <= 0:
        return None
    return (float((x*weights).sum()/total), float((y*weights).sum()/total))


def _finite(value):
    return None if not np.isfinite(value) else float(value)


def summarize(records, images, variances, refs, target, aperture, diagnostics=False):
    if len(records) < 3:
        return {"disposition": "insufficient_events", "accepted_events": len(records), "events": records}, {}
    stack = np.asarray(images)
    diff, noise = _combine(stack, np.asarray(variances))
    if not np.any(np.isfinite(diff) & aperture):
        return {"disposition": "insufficient_pixels", "accepted_events": len(records)}, {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        scene = np.nanmedian(refs, axis=0)
    peak = np.unravel_index(np.nanargmax(diff), diff.shape)
    center = _centroid(diff)
    snr = diff / noise
    # Bootstrap whole events, allowing the peak position to move as well.
    # This captures repeat-event scatter but omits PRF asymmetry/WCS errors.
    rng = np.random.default_rng(271828)
    centers = []
    for _ in range(128):
        idx = rng.integers(0, len(stack), len(stack))
        mean, _ = _combine(stack[idx], np.asarray(variances)[idx])
        c = _centroid(mean)
        if c is not None:
            centers.append(c)
    radius = None
    offset = None
    if center:
        offset = float(np.hypot(center[0] - target[0], center[1] - target[1]))
        if len(centers) >= 96:
            radius = float(np.percentile(np.hypot(np.array(centers)[:,0] - center[0],
                                                  np.array(centers)[:,1] - center[1]), 95))
    full_aperture = np.all(np.isfinite(diff[aperture])) and np.all(np.isfinite(scene[aperture]))
    depth = None
    if full_aperture and scene[aperture].sum() > 0:
        depth = float(1e6 * diff[aperture].sum() / scene[aperture].sum())
    disposition = "weak_or_unstable_difference"
    if (center and radius is not None and diff[peak] > 0 and snr[peak] >= 6 and radius <= 0.75):
        if 0 in peak or any(i == n-1 for i, n in zip(peak, diff.shape)):
            disposition = "stamp_edge_followup"
        elif offset <= 1.0:
            disposition = "target_pixel_vicinity"
        elif offset - radius >= 1.0:
            disposition = "offset_source_followup"
        else:
            disposition = "ambiguous_position"
    result = {"disposition": disposition, "accepted_events": len(records), "events": records,
              "peak_xy": [int(peak[1]), int(peak[0])], "centroid_xy": list(center) if center else None,
              "offset_pixels": offset, "event_bootstrap_radius95_pixels": radius,
              "peak_score": _finite(snr[peak]), "aperture_depth_ppm": depth,
              "complete_aperture": bool(full_aperture)}
    arrays = {"scene": scene, "difference": diff, "score": snr} if diagnostics else {}
    return result, arrays


def _lists(array):
    return [[_finite(x) for x in row] for row in array]


def followup(task, cache):
    """Analyze one TESS product and an independently supplied transit ephemeris."""
    return _followup(task, cache, False)


def variable_followup(task, cache):
    """Localize a frozen stellar period and minimum-light epoch in one TESS stamp."""
    period = float(task["period_days"])
    phased = {**task, "epoch_btjd": task["epoch_min_btjd"], "duration_days": .2*period}
    return _followup(phased, cache, True)


def _followup(task, cache, stellar):
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    name, period, epoch, duration = validate(task)
    path, uri, url, sha = _acquire(Path(cache), name, task.get("source_sha256"))
    with fits.open(path, memmap=True) as hdus:
        h, data = hdus[0].header, hdus["PIXELS"].data
        ht = hdus["PIXELS"].header
        if (h["TICID"] != task["tic_id"] or h["SECTOR"] != task["sector"]
                or ht.get("TIMESYS") != "TDB" or ht.get("TIMEUNIT") != "d"
                or ht["BJDREFI"] + ht.get("BJDREFF", 0) != 2457000):
            raise ValueError("Unexpected TPF identity, time units or reference frame")
        time = np.asarray(data["TIME"], dtype=float)
        good = np.isfinite(time) & (data["QUALITY"] == 0)
        t = time[good]
        if len(t) < 2 or not np.all(np.diff(t) > 0):
            raise ValueError("Not enough quality-zero, sorted TPF cadences")
        cadence = float(np.median(np.diff(t)))
        cube, err = data["FLUX"][good], data["FLUX_ERR"][good]
        mask = np.array(hdus["APERTURE"].data, dtype=int)
        ap = (mask & 2) != 0
        wcs = WCS(hdus["APERTURE"].header).celestial
        ra, dec = float(h["RA_OBJ"]), float(h["DEC_OBJ"])
        target = tuple(float(v) for v in wcs.world_to_pixel_values(ra, dec))
        if not np.all(np.isfinite(target)) or not ap.any():
            raise ValueError("Missing target sky position, celestial stamp WCS or science aperture")
        base = {"schema_version": 1, "method": METHOD, "source_sha256": sha, "archive_uri": uri,
                "archive_url": url, "tic_id": int(h["TICID"]), "sector": int(h["SECTOR"]),
                "camera": int(h["CAMERA"]), "ccd": int(h["CCD"]), "observation_id": name,
                "ra_deg": ra, "dec_deg": dec, "target_xy": list(target), "shape": list(ap.shape),
                "period_days": period, "epoch_btjd": epoch, "duration_days": duration,
                "time_system": "BJD_TDB - 2457000 (days)", "cadences_total": len(time),
                "cadences_accepted_quality_zero": len(t), "median_cadence_days": cadence,
                "first_btjd": float(t[0]), "last_btjd": float(t[-1]), "role": task.get("role", "followup")}
        tmag = _finite(h.get("TESSMAG", np.nan))
        frame_seconds = _finite(ht.get("INT_TIME", np.nan))
    if cube.ndim != 3 or cube.shape != err.shape or cube.shape[1:] != mask.shape:
        raise ValueError("Unexpected TPF image dimensions")
    window_scale, intervals = (period, VARIABLE_WINDOWS) if stellar else (duration, TRANSIT_WINDOWS)
    if stellar:
        base.update(method="tess-stellar-phase-difference-v1", epoch_min_btjd=epoch,
                    phase_intervals_min_and_references=[list(w) for w in intervals])
    nominal = event_images(t, cube, err, epoch, period, window_scale, cadence, intervals)
    result, arrays = summarize(*nominal, target, ap, diagnostics=True)
    scene = arrays.get("scene")
    result["tess_magnitude"] = tmag
    result["frame_integration_seconds"] = frame_seconds
    result["peak_scene_electrons_per_second"] = _finite(np.nanmax(scene)) if scene is not None else None
    result["detector_cautions"] = []
    if (center := result.get("centroid_xy")) is not None:
        sky = wcs.pixel_to_world(*center)
        reference = SkyCoord(ra * u.deg, dec * u.deg)
        result.update(centroid_ra_deg=float(sky.ra.deg), centroid_dec_deg=float(sky.dec.deg),
                      sky_separation_arcsec=float(reference.separation(sky).arcsec))
        dx, dy = center[0] - target[0], center[1] - target[1]
        result["offset_detector_xy_pixels"] = [dx, dy]
        # A very bright star can bleed along a CCD column. This heuristic
        # identifies a competing instrumental explanation for an offset.
        if tmag is not None and tmag <= 7 and abs(dx) <= .5 and abs(dy) >= 1:
            result["detector_cautions"].append("bright_target_column_aligned_offset")
    controls = []
    for phase in PHASE_OFFSETS:
        # Omit every cadence within a half-duration of the actual prediction so
        # control sidebands cannot accidentally measure the original target dip.
        dt = (t - epoch + period/2) % period - period/2
        keep = np.abs(dt) > duration/2
        off = event_images(t[keep], cube[keep], err[keep], epoch + phase * period,
                           period, window_scale, cadence, intervals)
        control, _ = summarize(*off, target, ap)
        controls.append({"phase_offset": phase, **{k:v for k,v in control.items() if k != "events"}})
    # Use the positive-phase control background, apply an artificial spatially
    # displaced box deficit in that window, and measure its recovered center.
    yy, xx = np.indices(ap.shape)
    x0 = target[0] + (2.3 if target[0] < (ap.shape[1]-1)/2 else -2.3)
    y0 = target[1] + (1.4 if target[1] < (ap.shape[0]-1)/2 else -1.4)
    injection = {"disposition": "insufficient_coverage"}
    if scene is not None and off[1] and 1 <= x0 < ap.shape[1]-1 and 1 <= y0 < ap.shape[0]-1:
        prf = np.exp(-((xx-x0)**2 + (yy-y0)**2)/(2*0.65**2))
        amplitude = 0.01 * max(0, float(np.nansum(scene[ap])))
        mock = [d + amplitude*prf/prf.sum() for d in off[1]]
        test, _ = summarize(off[0], mock, off[2], off[3], target, ap)
        recovered = test.get("centroid_xy")
        injection = {"injected_xy": [x0, y0], "amplitude_electrons_per_second": amplitude,
                     "recovered_xy": recovered, "disposition": test["disposition"],
                     "centroid_error_pixels": (float(np.hypot(recovered[0]-x0, recovered[1]-y0))
                                               if recovered else None)}
    return {**base, "disposition": result["disposition"], "localization": result,
            "controls": {"phase_shifts": controls, "off_target_injection": injection},
            "images": {k:_lists(v) for k,v in arrays.items()}, "science_aperture": ap.astype(int).tolist()}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("task", type=Path, help="Single task or JSON array of tasks")
    p.add_argument("--cache", type=Path, default=Path("runs/tess-pixels"))
    p.add_argument("--stellar", action="store_true", help="Measure bright-minus-faint stellar phase images")
    args = p.parse_args(argv)
    try:
        requested = json.loads(args.task.read_bytes())
        tasks = requested if isinstance(requested, list) else [requested]
        measure = variable_followup if args.stellar else followup
        rows = [measure(t, args.cache) for t in tasks]
        print(json.dumps(rows if isinstance(requested,list) else rows[0], indent=2, allow_nan=False))
    except (OSError, ValueError, KeyError, ImportError) as exc:
        p.exit(1, f"pixel follow-up: {exc}\n")


if __name__ == "__main__":
    main()
