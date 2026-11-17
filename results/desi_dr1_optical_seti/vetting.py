"""DESI DR1 candidate context tests. Wavelengths are vacuum Angstrom throughout."""

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SKY = tuple(float(s) for s in (HERE / "kpno_sky_vacuum.txt").read_text().splitlines()
            if s and not s.startswith("#"))
with (HERE / "rest_lines.csv").open() as f:
    REST = {r["name"]: float(r["vacuum_angstrom"]) for r in csv.DictReader(f)}
# Lamp lines tabulated by NIST in air. Invert the SDSS/Morton formula.
def air_to_vacuum(air):
    vac = air
    for _ in range(4):
        vac = air * (1 + 2.735182e-4 + 131.4182 / vac**2 + 2.76249e8 / vac**4)
    return vac


LAMP = {f"HgI_{air}": air_to_vacuum(air) for air in
        (4046.563, 4358.328, 5460.735, 5769.598, 5790.663)}
LAMP.update({name: REST[name] for name in ("nai_d2", "nai_d1")})
# Identified after the first unmasked stellar pass. Preserve the distinction
# from the pre-inspection rest_lines.csv atlas in every output record.
POST_INSPECTION = {"LiI_6707.775": air_to_vacuum(6707.775),
                   "LiI_6707.926": air_to_vacuum(6707.926)}
ANCHORS = ("oii_3726", "oii_3729", "oii_blend", "oiii_4959", "oiii_5007",
           "hbeta", "halpha", "lyalpha")
COMPANIONS = {name: wave for name, wave in REST.items()
              if name not in {"caii_k", "caii_h", "caii_8498", "caii_8542",
                              "caii_8662", "nai_d1", "nai_d2"}}
REVIEWS = ("stellar_template_and_transition_atlas", "detector_trace_and_cosmics",
           "sky_model_and_tellurics", "bright_neighbour_leakage",
           "nebula_and_emission_star", "galaxy_qso_and_imaging",
           "stellar_vs_terrestrial_frame")


def finite(value):
    return isinstance(value, (float, int)) and not isinstance(value, bool) and math.isfinite(value)


def positive(value):
    return finite(value) and value > 0


def reference_mask(bary_A, lsf_fwhm_A, topo_to_bary, star_to_bary, rv_error_kms,
                   include_post_inspection=True):
    """Return True where catalogue vetoes remove coverage, before any other tests."""
    import numpy as np

    wave, width = np.broadcast_arrays(np.asarray(bary_A, float), np.asarray(lsf_fwhm_A, float))
    if (not all(positive(x) for x in (topo_to_bary, star_to_bary))
            or not finite(rv_error_kms) or rv_error_kms < 0):
        raise ValueError("missing or invalid wavelength frame information")
    excluded = ~(np.isfinite(wave) & np.isfinite(width) & (width > 0))

    def apply(lines, center, radius):
        # The nearest catalogue centroid is enough to decide whether any window overlaps.
        ordered = np.sort(list(lines))
        right = np.searchsorted(ordered, center).clip(0, len(ordered) - 1)
        left = np.maximum(0, right - 1)
        return np.minimum(abs(center - ordered[left]), abs(center - ordered[right])) <= radius

    excluded |= apply((*SKY, *LAMP.values()), wave / topo_to_bary, width / topo_to_bary)
    rest = wave / star_to_bary
    lines = (*REST.values(), *POST_INSPECTION.values()) if include_post_inspection else REST.values()
    excluded |= apply(lines, rest, np.hypot(width / star_to_bary,
                                                      3 * rv_error_kms * rest / 299792.458))
    return excluded


def background_hypotheses(bary_wave, fwhm, other_peaks, coverage):
    """Project alternative redshifts, leaving one-line Ly-alpha and faint companions open.

    other_peaks are independently measured positive science features (barycentric A);
    coverage is the disjoint union of good-pixel intervals in that wavelength frame.
    """
    hypotheses = []
    for anchor in ANCHORS:
        # Equal components give a nominal unresolved [O II] centroid. The
        # individual anchors bracket the effect of an unknown doublet ratio.
        anchor_wave = 0.5 * (REST["oii_3726"] + REST["oii_3729"]) if anchor == "oii_blend" else REST[anchor]
        stretch = bary_wave / anchor_wave
        if stretch < 1:
            continue
        predictions = {}
        for name, rest in COMPANIONS.items():
            obs = stretch * rest
            if abs(obs - bary_wave) <= fwhm or not any(a <= obs <= b for a, b in coverage):
                continue
            predictions[name] = {
                "barycentric_A": round(obs, 3),
                "peaks_within_one_fwhm_A": [w for w in other_peaks
                                             if abs(w - obs) <= fwhm and abs(w - bary_wave) > fwhm],
            }
        hypotheses.append({"anchor": anchor, "z": stretch - 1, "companions": predictions})
    return hypotheses


def select_controls(fibers, is_sky, ra_deg, dec_deg, target_fiber, nsky=3):
    """Select unique sky fibres by angular distance and slit traces by fibre ID.

    Supply *all* fibres from one cframe camera/exposure, including flagged
    rows. A missing physical neighbour remains missing for adjudication.
    """
    import numpy as np

    fibers = np.asarray(fibers)
    is_sky = np.asarray(is_sky, dtype=bool)
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    if not (fibers.shape == is_sky.shape == ra.shape == dec.shape and fibers.ndim == 1):
        raise ValueError("malformed fibermap columns")
    if (np.unique(fibers).size != len(fibers) or target_fiber not in fibers
            or np.any(fibers // 500 != target_fiber // 500)):
        raise ValueError("expected one complete camera and a unique target fibre")
    i = np.flatnonzero(fibers == target_fiber)[0]
    hav = np.sin((dec - dec[i]) / 2)**2 + np.cos(dec) * np.cos(dec[i]) * np.sin((ra - ra[i]) / 2)**2
    candidates = np.flatnonzero(is_sky & (fibers != target_fiber) & np.isfinite(hav))
    sky = candidates[np.argsort(hav[candidates], kind="stable")[:nsky]]
    return {"sky_fibers": [int(fibers[j]) for j in sky],
            "sky_separation_arcsec": [float(np.rad2deg(2 * np.arcsin(np.sqrt(np.clip(hav[j], 0, 1)))) * 3600)
                                      for j in sky],
            "trace_fibers": [f for f in (target_fiber - 1, target_fiber + 1)
                             if f // 500 == target_fiber // 500]}


def forced_line(wave, flux, ivar, mask, line_profile):
    """Fit a local line and linear continuum on a cutout of one control trace.

    line_profile is the exposure/fibre's resolution-column density on this grid;
    supply flanking continuum bins on both sides and use the science flux scale.
    Returns integrated flux and its formal error or None for unusable data.
    """
    import numpy as np

    wave, flux, ivar, mask, p = map(np.asarray, (wave, flux, ivar, mask, line_profile))
    if not (wave.shape == flux.shape == ivar.shape == mask.shape == p.shape and wave.ndim == 1):
        raise ValueError("cutout columns disagree")
    if len(wave) < 7 or not np.all(np.isfinite(wave)) or np.any(np.diff(wave) <= 0):
        return None
    if not np.all(np.isfinite(p)) or max(p) <= 0:
        return None
    good = np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0) & (mask == 0)
    core = p >= 0.10 * max(p)
    below, above = np.flatnonzero(core)[[0,-1]]
    if (not np.all(good[core]) or good[:below].sum() < 2 or good[above+1:].sum() < 2):
        return None
    x = (wave - wave.mean()) / np.ptp(wave)
    design = np.column_stack((p[good], np.ones(good.sum()), x[good])) * np.sqrt(ivar[good])[:, None]
    y = flux[good] * np.sqrt(ivar[good])
    coef, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    if rank != 3:
        return None
    err = np.sqrt(np.linalg.inv(design.T @ design)[0, 0])
    return {"flux": float(coef[0]), "sigma": float(err)}


def verified_review(item, root):
    """Resolve an evidence file relative to the packet directory and check its digest."""
    if not isinstance(item, dict) or item.get("status") not in {"clear", "veto", "ambiguous"}:
        return None
    if not item.get("reason") or not item.get("file") or not isinstance(item.get("file"), str):
        return None
    path = (root / item["file"]).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file():
        return None
    if hashlib.sha256(path.read_bytes()).hexdigest() != item.get("sha256"):
        return None
    return item["status"]


def vet(packet, root=HERE):
    """Adjudicate one candidate packet; missing measurements always remain pending."""
    veto, pending, measurements = [], [], []
    expected = packet.get("expected_expids", [])
    epochs = packet.get("exposures", [])
    ids = [e.get("expid") for e in epochs]
    if (len(set(expected)) < 2 or len(set(expected)) != len(expected) or
            set(expected) != set(ids) or len(set(ids)) != len(ids) or None in ids):
        pending.append("complete_distinct_exposure_inventory")
    for e in epochs:
        label = str(e.get("expid", "?"))
        desc = lambda what: f"{label}/{what}"
        bary = e.get("barycentric_A")
        width = e.get("lsf_fwhm_A")
        h = e.get("topo_to_bary")
        doppler = e.get("star_to_bary")
        if not all(positive(v) for v in (bary, width, h, doppler)):
            pending.append(desc("wavelength_factors_or_width"))
            continue
        rest = bary / doppler
        topo = bary / h
        rv_error = e.get("rv_error_kms")
        if not finite(rv_error) or rv_error < 0:
            pending.append(desc("rv_uncertainty"))
            rv_error = 0.0
        tolerance = math.hypot(width / doppler, 3 * rv_error * rest / 299792.458)
        lines = [k for k, w in REST.items() if abs(w - rest) <= tolerance]
        post = [k for k, w in POST_INSPECTION.items() if abs(w - rest) <= tolerance]
        airglow = [w for w in SKY if abs(w - topo) <= width / h]
        lamps = [k for k, w in LAMP.items() if abs(w - topo) <= width / h]
        # Conservative wavelength exclusion. It does not identify the emitter.
        if lines:
            veto.append(desc("excluded_stellar_transition_window"))
        if post:
            veto.append(desc("post_inspection_stellar_transition_window"))
        if airglow or lamps:
            veto.append(desc("excluded_terrestrial_line_window"))
        sky = e.get("sky_controls", [])
        expected_sky = e.get("expected_sky_fibers", [])
        target = e.get("fiber")
        sky_ids = [x.get("fiber") for x in sky]
        if (not isinstance(target, int) or not 0 <= target < 5000 or len(set(expected_sky)) != 3
                or len(sky_ids) != 3 or set(sky_ids) != set(expected_sky) or target in sky_ids
                or any(not isinstance(f, int) or f // 500 != target // 500 for f in sky_ids)):
            pending.append(desc("three_nearest_distinct_sky_fibers"))
        for s in sky:
            amp, err = s.get("flux"), s.get("sigma")
            if not finite(amp) or not positive(err) or not positive(e.get("target_flux")):
                pending.append(desc("sky_flux_measurement"))
            elif abs(amp) >= 4 * err:
                veto.append(desc(f"shared_sky_residual_fiber_{s['fiber']}"))
            elif abs(amp) + 3 * err >= e["target_flux"]:
                pending.append(desc(f"insensitive_sky_fiber_{s['fiber']}"))
        if isinstance(target, int) and 0 <= target < 5000:
            neighbors = {f for f in (target-1, target+1) if f//500 == target//500}
            trace_ids = e.get("trace_fibers_measured", [])
            if set(trace_ids) != neighbors:
                pending.append(desc("adjacent_physical_traces"))
        coverage, peaks = e.get("good_bary_intervals_A"), e.get("other_peaks_bary_A")
        hypotheses = []
        if (isinstance(coverage, list) and coverage and isinstance(peaks, list)
                and all(isinstance(row, list) and len(row)==2 and all(positive(w) for w in row)
                        and row[0] < row[1] for row in coverage) and all(positive(w) for w in peaks)
                and any(a <= bary <= b for a,b in coverage)):
            hypotheses = background_hypotheses(bary, width, peaks, coverage)
        else:
            pending.append(desc("background_companion_coverage_or_peaks"))
        measurements.append({"expid": label, "stellar_rest_A": rest, "topocentric_A": topo,
                             "stellar_window_hits": lines, "post_inspection_stellar_hits": post,
                             "sky_window_hits_A": airglow,
                             "lamp_window_hits": lamps, "background_hypotheses": hypotheses})

    for check in REVIEWS + ("parent_inventory_persistence_and_lsf",):
        status = verified_review(packet.get("reviews", {}).get(check), Path(root))
        if status == "veto":
            veto.append(check)
        elif status != "clear":
            pending.append(check)
    return {"targetid": packet.get("targetid"), "status": "excluded" if veto else
            "pending" if pending else "unidentified_follow_up",
            "veto_reasons": sorted(set(veto)), "pending_checks": sorted(set(pending)),
            "measurements": measurements}


if __name__ == "__main__":
    packet_file = Path(sys.argv[1])
    print(json.dumps(vet(json.loads(packet_file.read_text()), packet_file.parent),
                     indent=2, allow_nan=False))
