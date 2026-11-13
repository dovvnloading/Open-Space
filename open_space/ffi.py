"""Search every interior pixel of a pair of TESS full-frame cutouts for recurring variability."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .ffi_archive import acquire
from .tess_archive import atomic_bytes
from .tess_eb_catalog import separation


METHOD = "tess-ffi-periodic-image-v1"
MIN_PERIOD, MAX_PERIOD = 0.18, 3.0
BORDER = 2
MIN_POWER, MIN_AMPLITUDE_SNR = 0.08, 8.0
MATCH_ARCSEC = 31.5


def _remove_slow(time, matrix):
    clean = np.empty_like(matrix)
    for idx in np.split(np.arange(len(time)), np.flatnonzero(np.diff(time) > 0.3) + 1):
        tx, values = time[idx], matrix[idx]
        anchors = np.linspace(tx[0], tx[-1], max(2, int(np.ceil(np.ptp(tx)/0.25))+1))
        baseline = np.stack([np.median(values[abs(tx-u) <= 0.9], axis=0) for u in anchors])
        left = np.clip(np.searchsorted(anchors, tx, side="right") - 1, 0, len(anchors)-2)
        w = ((tx-anchors[left])/np.maximum(anchors[left+1]-anchors[left], 1e-12))[:,None]
        clean[idx] = values - baseline[left]*(1-w) - baseline[left+1]*w
    return clean - clean.mean(axis=0)


def condition(time, images, errors, excluded=None):
    """Remove diffuse backgrounds, small translations and changes in common illumination."""
    x_count = images.shape[2]
    matrix = images.reshape(len(time), -1)
    error = errors.reshape(len(time), -1)
    available = (np.isfinite(matrix) & np.isfinite(error) & (error > 0)).mean(axis=0) >= 0.99
    if available.sum() < 60:
        raise ValueError("Too few stable pixels for field-background estimation")
    reference = np.nanmedian(matrix, axis=0)
    reference = np.where(np.isfinite(reference), reference, 0.)
    matrix = np.where(np.isfinite(matrix), matrix, reference)
    matrix = matrix - reference
    profile = reference.reshape(images.shape[1:])
    gy, gx = np.gradient(profile)
    yy, xx = np.indices(profile.shape, dtype=float)
    spatial = np.column_stack([np.ones(reference.size), xx.ravel()-xx.mean(), yy.ravel()-yy.mean(),
                               gx.ravel(), gy.ravel(), reference - np.median(reference)])
    fit_pixels = available & (reference <= np.quantile(reference[available], 0.80))
    if excluded is not None:
        x, y = excluded
        fit_pixels &= ((abs(xx-x) > 3) | (abs(yy-y) > 3)).ravel()
    if fit_pixels.sum() < 20:
        raise ValueError("Too few background pixels outside the candidate footprint")
    design = spatial[fit_pixels]
    coeff = matrix[:,fit_pixels] @ np.linalg.pinv(design).T
    # Limit the leverage of faint variable sources and moving objects in a frame's nuisance fit.
    residuals = matrix[:,fit_pixels] - coeff @ design.T
    scatter = np.maximum(1e-5, 1.4826*np.median(abs(residuals-np.median(residuals,axis=1)[:,None]),axis=1))
    for j in range(len(time)):
        keep = abs(residuals[j]) < 3.5*scatter[j]
        if keep.sum() >= 20:
            coeff[j] = np.linalg.lstsq(design[keep], matrix[j,fit_pixels][keep], rcond=None)[0]
    corrected = _remove_slow(time, matrix-coeff @ spatial.T)
    # Match the search's high pass when comparing against detector-motion timescales.
    regressors = _remove_slow(time, coeff)
    good = np.ones(profile.shape, dtype=bool)
    good[:BORDER] = good[-BORDER:] = good[:,:BORDER] = good[:,-BORDER:] = False
    return corrected, good.ravel() & available, regressors, x_count


def search_frequencies(time, series):
    """Evaluate a floating-mean sine and cosine at each trial frequency for all pixels."""
    time = time - time[0]
    frequencies = np.linspace(1/MAX_PERIOD, 1/MIN_PERIOD,
                              int(np.ceil((1/MIN_PERIOD-1/MAX_PERIOD)*np.ptp(time)*8))+1)
    variance = np.maximum(1e-20, np.sum(series**2, axis=0))
    best = np.zeros(series.shape[1])
    best_frequency = np.full(series.shape[1], frequencies[0])
    amplitudes = np.zeros_like(best)
    cols = np.arange(series.shape[1])
    for group in np.array_split(frequencies, max(1, int(np.ceil(len(frequencies)/128)))):
        angle = 2*np.pi*group[:,None]*time
        sine, cosine = np.sin(angle), np.cos(angle)
        sine -= sine.mean(axis=1,keepdims=True)
        cosine -= cosine.mean(axis=1,keepdims=True)
        ss, cc, sc = (sine*sine).sum(axis=1), (cosine*cosine).sum(axis=1), (sine*cosine).sum(axis=1)
        det = (ss*cc-sc**2)[:,None]
        sy, cy = sine @ series, cosine @ series
        bs = (cc[:,None]*sy-sc[:,None]*cy)/det
        bc = (ss[:,None]*cy-sc[:,None]*sy)/det
        power = (bs*sy+bc*cy)/variance
        j = power.argmax(axis=0)
        update = power[j,cols] > best
        best[update] = power[j,cols][update]
        best_frequency[update] = group[j[update]]
        amplitudes[update] = np.hypot(bs[j,cols],bc[j,cols])[update]
    return best_frequency, best, amplitudes, len(frequencies)


def _wave(time, frequency):
    angle = 2*np.pi*frequency*(time-time[0])
    return np.column_stack([np.cos(angle), np.sin(angle), np.ones(len(time))])


def detect(cube):
    time, images, errors, wcs = cube
    corrected, use, nuisance, width = condition(time, images, errors)
    frequency, power, amplitude, trials = search_frequencies(time, corrected)
    noise = np.maximum(1e-8, (1.4826/np.sqrt(2))*np.median(abs(np.diff(corrected,axis=0)),axis=0))
    error_floor = np.nanmedian(errors,axis=0).ravel()
    snr = amplitude*np.sqrt(len(time)/2)/np.maximum(noise,error_floor)
    pass_pixel = use & (power >= MIN_POWER) & (snr >= MIN_AMPLITUDE_SNR)
    # A variable's flux is shared between neighboring pixels. Find the brightest
    # spatial peak meeting the detection rule and suppress its immediate footprint.
    ranks = np.argsort(np.where(pass_pixel, amplitude*np.sqrt(power), 0.))[::-1]
    yy, xx = np.indices(images.shape[1:])
    sources = []
    for k in ranks:
        if not pass_pixel[k]:
            break
        y, x = divmod(int(k), width)
        if any(max(abs(s["peak_x"]-x),abs(s["peak_y"]-y)) <= 2 for s in sources):
            continue
        design = _wave(time,frequency[k])
        independent, _, independent_nuisance, _ = condition(time, images, errors, excluded=(x,y))
        fit = np.linalg.lstsq(design, independent, rcond=None)[0][:2]
        direction = fit[:,k] / max(1e-12, np.linalg.norm(fit[:,k]))
        projected = (direction @ fit).reshape(images.shape[1:])
        valid = np.isfinite(projected)
        patch = valid & (abs(xx-x)<=2) & (abs(yy-y)<=2)
        inner = (abs(xx-x)<=1) & (abs(yy-y)<=1)
        # Keep flux attributed to the central PSF by clipping the local background.
        surround = valid & ~inner & (abs(xx-x)<=3) & (abs(yy-y)<=3)
        local = projected - np.median(projected[surround])
        weights = np.maximum(local,0)*patch
        total = weights.sum()
        if total <= 0:
            continue
        cx,cy = (weights*xx).sum()/total,(weights*yy).sum()/total
        ra,dec = wcs.pixel_to_world_values(cx,cy)
        dipole = np.maximum(-local[patch],0).sum()/total
        concentration = weights[inner].sum()/total
        yseries = independent[:,k]
        denom = np.linalg.norm(independent_nuisance,axis=0)*np.linalg.norm(yseries)
        corr = np.max(abs(yseries @ independent_nuisance)/np.maximum(denom,1e-15))
        rechecked_power = np.sum((design[:,:2] @ fit[:,k] -
                                  np.mean(design[:,:2] @ fit[:,k]))**2) / max(1e-20,yseries @ yseries)

        flags = []
        if rechecked_power < MIN_POWER:
            flags.append("signal_lost_on_independent_nuisance_fit")
        if dipole > 0.55:
            flags.append("bipolar_residual")
        if concentration < 0.40 or np.hypot(cx-x,cy-y) > 1.5:
            flags.append("diffuse_spatial_response")
        if corr > 0.55:
            flags.append("field_nuisance_correlation")
        sources.append({"peak_x": x, "peak_y": y, "centroid_x": float(cx), "centroid_y": float(cy),
                        "ra_deg": float(ra), "dec_deg": float(dec), "period_days": float(1/frequency[k]),
                        "variance_fraction": float(power[k]), "semiamplitude_e_per_s": float(amplitude[k]),
                        "amplitude_snr_proxy": float(snr[k]), "negative_to_positive_flux": float(dipole),
                        "independent_nuisance_fit_variance_fraction": float(rechecked_power),
                        "core_flux_fraction": float(concentration), "max_field_regressor_correlation": float(corr),
                        "artifact_flags": flags})
    meta = {"eligible_pixels": int(use.sum()), "above_threshold_pixels": int(pass_pixel.sum()),
            "spatial_peaks": len(sources), "frequency_trials": trials, "valid_cadences": len(time)}
    return sources, meta


def reconcile(first, later):
    output = []
    matches = {}
    pairs = sorted((separation(s["ra_deg"], s["dec_deg"], q["ra_deg"], q["dec_deg"]),i,j)
                   for i,s in enumerate(first) for j,q in enumerate(later))
    claimed = set()
    for distance, i, j in pairs:
        if distance <= MATCH_ARCSEC and i not in matches and j not in claimed:
            matches[i] = (distance,j)
            claimed.add(j)
    for i, s in enumerate(first):
        distance, j = matches.get(i, (None, None))
        q = None if j is None else later[j]
        frac = None if q is None else abs(q["period_days"]/s["period_days"]-1)
        if s["artifact_flags"] or (q is not None and q["artifact_flags"]):
            status = "image_artifact_veto"
        elif q is None:
            status = "single_sector_image_peak"
        elif frac > 0.025:
            status = "frequency_inconsistent"
        else:
            status = "recurrent_image_variable"
        output.append({"source_id": f"ICRS-{s['ra_deg']:.5f}{s['dec_deg']:+.5f}",
                       "ra_deg": s["ra_deg"], "dec_deg": s["dec_deg"], "discovery": s, "validation": q,
                       "match_separation_arcsec": distance, "period_fractional_difference": frac,
                       "disposition": status, "candidate": status == "recurrent_image_variable"})
    return output, [q for j,q in enumerate(later) if j not in claimed]


def run_field(field, cache_dir):
    products = field["products"]
    if len(products) != 2 or products[0]["sector"] >= products[1]["sector"]:
        raise ValueError("Plan exactly two increasing sectors per field")
    searches, metadata = [], []
    for product in products:
        cube, meta = acquire(field, product, cache_dir)
        detections, counts = detect(cube)
        searches.append(detections)
        metadata.append({**meta, **counts})
    sources, unmatched = reconcile(*searches)
    return {"field_id": field["field_id"], "role": field["role"], "method": METHOD,
            "ra_deg": field["ra_deg"], "dec_deg": field["dec_deg"], "size_pixels": field["size_pixels"],
            "products": metadata, "sources": sources, "unmatched_validation_peaks": unmatched}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=Path(__file__).resolve().parents[1]/"data/ffi-fields.json")
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.plan.read_bytes()
    fields = [run_field(item,args.cache) for item in json.loads(raw)["fields"]]
    result = {"schema_version": 1, "method": METHOD, "plan_sha256": sha256(raw).hexdigest(), "fields": fields}
    atomic_bytes(args.output, (json.dumps(result, indent=2, allow_nan=False)+"\n").encode())


if __name__ == "__main__":
    main()
