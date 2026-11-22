"""Blind HETDEX 1D search. Input wavelengths are observed air Angstroms.

Input flux densities and 1-sigma errors are in physical cgs per Angstrom by
default. Pass density_unit=1e-17 for the public FITS numeric convention.
All reported fluxes/flux densities are in physical cgs. No catalog metadata
enter detection. See METHODS.md for interpretation of the screening flags.
"""
from __future__ import annotations

import math
import warnings

import numpy as np
from scipy.ndimage import convolve1d, maximum_filter1d, uniform_filter1d


MIN_SNR = 7.0
SCREEN_SNR = 5.5
BANK = (2.0, 3.7, 6.5)
WIDTHS = (1.8, 2.4, 3.5, 5.0, 7.5, 10.0)
C = 299792.458
# Vacuum wavelengths; doublet components have separate names.
UV = {"lya": 1215.6701, "nv1": 1238.821, "nv2": 1242.804,
      "civ1": 1548.204, "civ2": 1550.781, "heii": 1640.42,
      "oiii_uv1": 1660.809, "oiii_uv2": 1666.150,
      "ciii1": 1906.683, "ciii2": 1908.734}
OPT = {"oii1": 3727.092, "oii2": 3729.875, "neiii": 3869.86,
       "hdelta": 4102.892, "hgamma": 4341.684, "hbeta": 4862.683,
       "oiii4959": 4960.295, "oiii5007": 5008.240, "halpha": 6564.61}
# Native air wavelengths: flags are warnings about coincidence, not subtractions.
SKY_AIR = np.array([4046.56, 4358.33, 5197.90, 5200.26, 5460.74, 5577.34])
# Published HETDEX intervals with bad wavelength-dependent calibration.
BAD_AIR = ((3534, 3556), (5194, 5197), (5200, 5205), (5456, 5466))


def air_to_vac(w):
    """Standard dry-air optical index (Morton); native grid remains unshifted."""
    s2 = (1e4 / np.asarray(w)) ** 2
    return np.asarray(w) * (1 + 8.34254e-5 + 0.02406147 / (130 - s2)
                            + 0.00015998 / (38.9 - s2))


def vac_to_air(v):
    v = float(v)
    air = v / 1.00028
    return float(air * v / air_to_vac(air))


def _weighted_cont(y, e, minimum=8):
    """Robust centered sideband and its error, with a block-dispersion floor."""
    ok = np.isfinite(y) & np.isfinite(e) & (e > 0)
    for low, high in BAD_AIR:
        ok[:, (w >= low) & (w <= high)] = False
    if np.count_nonzero(ok) < minimum:
        return None
    # Contiguous 8-pixel blocks retain local spectral correlations. Compute
    # the formal floor from exactly the same weights used in the mean.
    center = np.median(y[ok])
    scatter = 1.4826 * np.median(np.abs(y[ok] - center))
    ok &= np.abs(y - center) <= 4.5 * np.maximum(e, scatter)
    if np.count_nonzero(ok) < minimum:
        return None
    weights = np.zeros_like(e)
    weights[ok] = 1 / e[ok] ** 2
    sums = np.add.reduceat(weights, np.arange(0, len(y), 8))
    products = np.multiply(weights, y, out=np.zeros_like(y), where=ok)
    values = np.add.reduceat(products, np.arange(0, len(y), 8))
    mean = np.sum(values) / np.sum(sums)
    formal = 1 / np.sqrt(np.sum(sums))
    keep = sums > 0
    if np.count_nonzero(keep) >= 3:
        # Cluster-robust variance for the weighted mean, treating wavelength
        # blocks as independent. Assumes neither extraction nor calibration
        # introduces appreciable correlations beyond this scale.
        n = np.count_nonzero(keep)
        spread = np.sqrt(n / (n - 1) * np.sum((values[keep] - mean * sums[keep])**2))
        formal = max(formal, spread / np.sum(sums))
    return float(mean), float(formal), int(np.count_nonzero(ok))


def _grid(wave):
    w = np.asarray(wave, dtype=np.float64)
    if w.ndim != 1 or len(w) < 64 or not np.all(np.isfinite(w)):
        raise ValueError("Expected the finite, one-dimensional HETDEX air wavelength grid")
    dx = float(np.median(np.diff(w)))
    if dx <= 0 or np.max(np.abs(np.diff(w) - dx)) > dx * 0.005:
        raise ValueError("Expected a uniformly increasing, native HETDEX wavelength grid")
    return w, dx


def _inputs(w, flux, err, density_unit):
    y = np.atleast_2d(np.asarray(flux, dtype=np.float64)) * (density_unit * 1e17)
    e = np.atleast_2d(np.asarray(err, dtype=np.float64)) * (density_unit * 1e17)
    if y.shape != e.shape or y.shape[1] != len(w):
        raise ValueError("Matched flux, error, wavelength dimensions are required")
    ok = np.isfinite(y) & np.isfinite(e) & (e > 0)
    # Reject catalog masked bins with invalid variance, leaving real zero
    # flux with a positive variance as a valid measurement.
    return np.where(ok, y, 0), np.where(ok, e, np.inf), ok


def _screen(y, e, ok, dx):
    weights = 1 / e**2
    span = 2 * int(round(65 / dx)) + 1
    smooth_w = uniform_filter1d(weights, span, axis=1, mode="constant")
    smooth_y = uniform_filter1d(y * weights, span, axis=1, mode="constant")
    base = np.divide(smooth_y, smooth_w, out=np.zeros_like(y), where=smooth_w > 0)
    residual = (y - base) * weights
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        # Nominal errors remain the floor. Estimate marginal scatter from the
        # central same-spectrum channels after removing the local continuum.
        signed = np.where(ok, (y - base) / e, np.nan)
        mid = np.nanmedian(signed, axis=1, keepdims=True)
        noise = np.nanmedian(np.abs(signed - mid), axis=1) * 1.4826
    scale = np.where(np.isfinite(noise), np.maximum(noise, 1.0), 1)
    best = np.zeros_like(y)
    for sigma in BANK:
        r = int(math.ceil(3.5 * sigma / dx))
        kernel = np.exp(-0.5 * (np.arange(-r, r + 1) * dx / sigma)**2)
        num = convolve1d(residual, kernel, axis=1, mode="constant")
        den = convolve1d(weights, kernel**2, axis=1, mode="constant")
        sn = np.divide(num, np.sqrt(den), out=np.zeros_like(y), where=den > 0)
        np.maximum(best, sn, out=best)
    best /= scale[:, None]
    best[~ok] = 0
    peaks = (best >= SCREEN_SNR) & (maximum_filter1d(best, 5, axis=1) == best)
    return best, np.argwhere(peaks), scale


def _slice(w, left, right):
    return slice(int(np.searchsorted(w, left)), int(np.searchsorted(w, right, side="right")))


def _local_side(w, y, e, center, width):
    gap = max(12.0, 3.8 * width)
    side = _slice(w, center - gap - 48, center + gap + 48)
    if side.start >= side.stop:
        return None
    es = np.where(np.abs(w[side] - center) >= gap, e[side], np.inf)
    return _weighted_cont(y[side], es)


def _profile(w, y, e, center, sigma, continuum, dx):
    """Integrated Gaussian area, with sideband-error propagation."""
    if continuum is None:
        return None
    win = _slice(w, center - 4.1*sigma, center + 4.1*sigma)
    off = w[win] - center
    if len(off) < 3:
        return None
    g = np.exp(-0.5 * (off / sigma)**2) / (math.sqrt(2*math.pi) * sigma)
    wt = 1 / e[win]**2
    # All profile area, including outside the recorded band, belongs in
    # the coverage denominator. Missing central samples are never imputed.
    coverage = np.sum(g[wt > 0]) * dx
    if coverage < 0.82 or np.count_nonzero(wt[g >= g.max()/2]) < 2:
        return None
    den = float(np.dot(wt, g*g))
    if den <= 0:
        return None
    numerator = np.dot(wt*g, y[win] - continuum[0])
    error = math.hypot(1/math.sqrt(den), continuum[1]*np.dot(wt, g)/den)
    return float(numerator / den), float(error)


def _fit(w, y, e, dx, center):
    # Use one common sideband for model width comparison so absorption or
    # continuum differences cannot be maximized over the width bank.
    side = _local_side(w, y, e, center, max(WIDTHS))
    if side is None:
        return None
    best = None
    for sigma in WIDTHS:
        pair = _profile(w, y, e, center, sigma, side, dx)
        if pair is None:
            continue
        area, error = pair
        sn = area/error
        if best is None or sn > best[0]:
            best = (sn, area, error, sigma)
    if best is None or best[0] < MIN_SNR:
        return None
    sn, area, error, sigma = best
    win = _slice(w, center-2.5*sigma, center+2.5*sigma)
    pos = np.maximum(0, y[win] - side[0]) * np.isfinite(e[win])
    one = bool(pos.sum() > 0 and pos.max()/pos.sum() > 0.68)
    vac = float(air_to_vac(center))
    zl = vac / UV["lya"] - 1
    zh = vac / UV["heii"] - 1
    # Red-side continuum only (observed-frame interval); no assignment of
    # infinite EW to a non-detection. Require adequate wavelength coverage.
    red = _slice(w, center + max(25, 4*sigma), center + 245)
    continuum = _weighted_cont(y[red], e[red], minimum=32)
    ew = ew_lower = cont = cont_err = None
    if continuum is not None:
        cont, cont_err, _ = continuum
        if cont > 0:
            ew = area/(cont*(1+zl))
        # Two simultaneous one-sided intervals (Bonferroni joint >=95%
        # for an a-priori line under the Gaussian extraction noise model).
        # Physical continuum clamp avoids infinite or negative limits.
        upper = max(cont, 0) + 1.96*cont_err
        ew_lower = max(area - 1.96*error, 0)/(upper*(1+zl))
    return {"wave": float(center), "wave_vac": vac, "flux": area*1e-17,
            "flux_err": error*1e-17, "snr": sn, "sigma_A": sigma,
            "sky_near": bool(np.min(np.abs(SKY_AIR-center)) <= max(6, 2*sigma)),
            "single_pixel": one, "z_lya": zl, "z_heii": zh if 1.9 <= zh <= 2.4 else None,
            "ew0": ew, "ew0_lower": ew_lower, "ew_extreme": bool(ew_lower is not None and ew_lower > 240),
            "continuum": None if cont is None else cont*1e-17,
            "continuum_err": None if cont_err is None else cont_err*1e-17}


def search_batch(wave, flux, err, *, density_unit=1.0):
    """Yield list of (input row offset, peak dict); process ~128–512 rows/call."""
    w, dx = _grid(wave)
    y, e, ok = _inputs(w, flux, err, density_unit)
    score, peaks, scale = _screen(y, e, ok, dx)
    result = []
    for row, pix in peaks:
        center = w[pix]
        if 0 < pix < len(w)-1:
            l, m, r = score[row, pix-1:pix+2]
            curvature = l - 2*m + r
            if curvature < -1e-8:
                center += dx * float(np.clip((l-r)/(2*curvature), -0.5, 0.5))
        peak = _fit(w, y[row], e[row]*scale[row], dx, center)
        if peak is not None:
            peak["noise_scale"] = float(scale[row])
            peak["lya_pre_snr"] = None
            if peak["z_heii"] is not None:
                predicted = vac_to_air(UV["lya"]*(1+peak["z_heii"]))
                near = _slice(w, predicted*(1-300/C)-dx/2, predicted*(1+600/C)+dx/2)
                if near.stop > near.start and np.any(ok[row, near]):
                    peak["lya_pre_snr"] = float(np.max(score[row, near]))
            result.append((int(row), peak))
    return result


def search(wave, flux, err, *, density_unit=1.0):
    return [p for _, p in search_batch(wave, flux, err, density_unit=density_unit)]


def _force(w, y, e, dx, vac, sigma):
    center = vac_to_air(vac)
    if center < w[0] or center > w[-1]:
        return {"status": "outside_band"}
    side = _local_side(w, y, e, center, sigma)
    pair = _profile(w, y, e, center, sigma, side, dx)
    if pair is None:
        return {"status": "insufficient_coverage"}
    f, se = pair
    return {"status": "measured", "air_wave": center, "flux": f*1e-17,
            "flux_err": se*1e-17, "limit3": (max(f, 0)+3*se)*1e-17, "snr": f/se}


def detail(wave, flux, err, peak, *, density_unit=1.0):
    """Forced component fluxes/3σ limits at the observed peak redshift.

    Components are fit individually; blended components are correlated,
    so independent component limits must not be added in quadrature.
    Offsets are fixed except for the reported Ly-alpha velocity search.
    """
    w, dx = _grid(wave)
    y, e, _ = _inputs(w, flux, err, density_unit)
    y, e = y[0], e[0]*peak["noise_scale"]
    sigma = peak["sigma_A"]
    vac = peak["wave_vac"]
    result = {}
    if peak["z_heii"] is not None:
        z1 = 1+peak["z_heii"]
        family = {k: _force(w, y, e, dx, rest*z1, sigma) for k, rest in UV.items()}
        trials = []
        for vel in (-300, 0, 300, 600):
            line = _force(w, y, e, dx, UV["lya"]*z1*(1+vel/C), sigma)
            if line["status"] == "measured":
                trials.append(dict(line, offset_km_s=vel))
        result["heii_family"] = family
        result["lya_velocity_best"] = max(trials, key=lambda p:p["snr"], default=None)
    zl1 = 1+peak["z_lya"]
    result["lya_family"] = {k: _force(w, y, e, dx, rest*zl1, sigma)
                            for k, rest in UV.items() if k != "lya"}
    # For each plausible low-z single-line identity explicitly test other
    # available optical transitions at the same redshift.
    low = {}
    for primary, rest in {"oii": 3728.48, "neiii": OPT["neiii"],
                          "hdelta": OPT["hdelta"], "hgamma": OPT["hgamma"],
                          "hbeta": OPT["hbeta"], "oiii4959": OPT["oiii4959"],
                          "oiii5007": OPT["oiii5007"], "halpha": OPT["halpha"]}.items():
        z1 = vac/rest
        if z1 < 0.999:
            continue
        companions = {k: _force(w, y, e, dx, r*z1, sigma) for k, r in OPT.items()
                      if k != primary and not (primary == "oii" and k.startswith("oii"))}
        low[primary] = {"z": z1-1, "companions": companions}
    result["low_z"] = low
    result["unidentified"] = False  # Ly-alpha alone fits every wavelength in this band.
    return result


def shuffle_fibres(flux_matrix, err_matrix, seed):
    """Swap row identities independently at each wavelength. Pair flux/error.

    For sky-residual structure to persist, group by shot/extraction class
    wherever possible. This does not preserve per-row spectral covariance.
    """
    y, e = np.asarray(flux_matrix), np.asarray(err_matrix)
    if y.ndim != 2 or y.shape != e.shape or len(y) < 2:
        raise ValueError("Need at least two same-grid fibres with paired errors")
    index = np.random.default_rng(seed).permuted(np.broadcast_to(np.arange(len(y))[:, None], y.shape), axis=0)
    return np.take_along_axis(y, index, axis=0), np.take_along_axis(e, index, axis=0)


def _inject(w, y, center, area, width):
    y += area * np.exp(-0.5*((w-center)/width)**2)/(math.sqrt(2*math.pi)*width)


def simulate(wave, flux, err, seed, n_inject=3, *, density_unit=1.0):
    """One copy per trial of a real spectrum. Preserve original noise and masks.

    Cycle HeII+Lyα, Lyα with controlled continuum, HeII+Lyα+CIV contaminants;
    use 5, 7, 10, 15 and 25 nominal sigma and varied observed widths.
    Lyα trials reset the measured local continuum to a nominal rest EW of
    120, 240, 360, 600, or 1000 Å by adding a constant flux density to the
    real noise and other lines. This is conditional on the baseline estimate.
    Multi-line recovery includes a
    forced Lyα S/N>=5 trial and an explicit metal veto at CIV S/N>=3.
    """
    w, dx = _grid(wave)
    y, e, _ = _inputs(w, flux, err, density_unit)
    y, e = y[0], e[0]
    rng = np.random.default_rng(seed)
    before = search(w, y, e, density_unit=1e-17)
    result = []
    for j in range(n_inject):
        family = ("heii_lya", "lya", "heii_lya_civ")[(int(seed)+j) % 3]
        lo, hi = ((max(w[0]+12, vac_to_air(UV["heii"]*2.90)), min(w[-1]-12, vac_to_air(UV["heii"]*3.40)))
                  if family != "lya" else (w[0]+60, w[-1]-65))
        if hi <= lo:
            continue
        center = rng.uniform(lo, hi)
        width = float(rng.choice((2.0, 3.5, 6.0)))
        target = float(rng.choice((5, 7, 10, 15, 25)))
        side = _local_side(w, y, e, center, width)
        measure = _profile(w, y, e, center, width, side, dx)
        row = {"family": family, "wave": center, "sigma_A": width, "target_snr": target,
               "was_blended": any(abs(p["wave"]-center) <= max(4, width) for p in before),
               "input_flux": None, "recovered": False, "selected_ew": False,
               "lya_recovered": None, "civ_detected": None, "status": "unusable"}
        if not measure:
            result.append(row)
            continue
        area = target*measure[1]
        inj = y.copy()
        _inject(w, inj, center, area, width)
        red = _slice(w, center+max(25, 4*width), center+245)
        base_c = _weighted_cont(y[red], e[red], minimum=32)
        r1 = float(air_to_vac(center)) / (UV["lya"] if family == "lya" else UV["heii"])
        row["input_ew0"] = None
        if family == "lya" and base_c is not None:
            row["input_ew0"] = float(rng.choice((120, 240, 360, 600, 1000)))
            inj += area/(row["input_ew0"]*r1) - base_c[0]
        if family != "lya":
            _inject(w, inj, vac_to_air(UV["lya"]*r1), 3*area, width)
            if family == "heii_lya_civ":
                _inject(w, inj, vac_to_air(UV["civ1"]*r1), area*0.67, width)
                _inject(w, inj, vac_to_air(UV["civ2"]*r1), area*0.33, width)
        matching = [p for p in search(w, inj, e, density_unit=1e-17)
                    if abs(p["wave"]-center) <= max(4, width)]
        row.update(status="injected", input_flux=area*1e-17, recovered=bool(matching))
        if matching:
            p = max(matching, key=lambda q:q["snr"])
            row["selected_ew"] = p["ew_extreme"]
            if family != "lya" and p["z_heii"] is not None:
                z1 = 1+p["z_heii"]
                errors = e*p["noise_scale"]
                lya = _force(w, inj, errors, dx, UV["lya"]*z1, p["sigma_A"])
                civ = [_force(w, inj, errors, dx, UV[k]*z1, p["sigma_A"]) for k in ("civ1", "civ2")]
                row["lya_recovered"] = lya.get("snr", 0) >= 5
                row["civ_detected"] = any(line.get("snr", 0) >= 3 for line in civ)
        result.append(row)
    return {"injections": result}
