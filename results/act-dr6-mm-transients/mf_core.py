"""Map-space weighted point-source matched-filter **kernel**, not a validated survey.

Input maps must be on the *same pixel grid*, in identical native units.
The reference must exclude the tested observation. Pixel ivars model WHITE noise
only. Correlated scan noise, beams, false-discovery trials, completeness and
catalog vetoes must be validated on actual data before astrophysical claims.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.ndimage import maximum_filter
from scipy.signal import fftconvolve


def gaussian_peak_psf(sigma_pix: float, radius_sigma: float = 4.0) -> np.ndarray:
    """Synthetic-test PSF only; do NOT substitute for measured ACT beam."""
    if not (math.isfinite(sigma_pix) and sigma_pix > 0):
        raise ValueError("sigma_pix must be positive and finite")
    r = math.ceil(radius_sigma * sigma_pix)
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return np.exp(-(x*x + y*y) / (2 * sigma_pix*sigma_pix))


def inject(arr: np.ndarray, psf: np.ndarray, amplitude: float, y: int, x: int) -> np.ndarray:
    """Inject a synthetic point source with known peak amplitude into a map."""
    a = np.asarray(arr, dtype=float).copy()
    p = np.asarray(psf, dtype=float)
    if p.ndim != 2 or min(p.shape) % 2 != 1 or max(p.shape) % 2 != 1:
        raise ValueError("PSF must have odd dimensions")
    hy, hx = p.shape[0] // 2, p.shape[1] // 2
    if y < hy or x < hx or y + hy >= a.shape[0] or x + hx >= a.shape[1]:
        raise ValueError("Injection overlaps map edge")
    a[y-hy:y+hy+1, x-hx:x+hx+1] += amplitude * p
    return a


def matched_filter(
    observed: np.ndarray,
    reference: np.ndarray,
    observed_ivar: np.ndarray,
    reference_ivar: np.ndarray,
    psf_peak_one: np.ndarray,
    *,
    independent_reference: bool = False,
    minimum_support: float = 0.95,
) -> dict[str, np.ndarray]:
    """Weighted least squares source amplitude on a map-minus-reference.

    rho(x) = sum_p PSF(p-x) * residual(p) * weight(p)
    kappa(x) = sum_p PSF(p-x)^2 * weight(p)
    weight = 1/(1/ivar_obs + 1/ivar_reference), for INDEPENDENT maps.
    Returns amplitude=rho/kappa in native map units and its WHITE-noise SNR.
    Invalid and map-edge positions are NaN. No correlated-noise model is implied.
    """
    if not independent_reference:
        raise ValueError("Reference contains or may contain observation: leave-one-out coadd required")
    if not (0 < minimum_support <= 1):
        raise ValueError("minimum_support must lie in (0,1]")
    arrays = [np.asarray(a, dtype=float) for a in (observed, reference, observed_ivar, reference_ivar)]
    o, r, io, ir = arrays
    if any(a.ndim != 2 or a.shape != o.shape for a in arrays):
        raise ValueError("All maps/ivars must have the identical 2D pixel grid")
    p = np.asarray(psf_peak_one, dtype=float)
    if p.ndim != 2 or any(d % 2 == 0 for d in p.shape) or not np.all(np.isfinite(p)):
        raise ValueError("PSF must be finite, odd-sized, two-dimensional")
    if not (np.isclose(np.max(p), 1.0, atol=1e-9) and np.all(p >= 0) and p.sum() > 0):
        raise ValueError("PSF must be nonnegative and peak-normalized to 1")
    valid = np.isfinite(o) & np.isfinite(r) & np.isfinite(io) & np.isfinite(ir) & (io > 0) & (ir > 0)
    # The where prevents NaNs and infs from contaminating the FFT convolution.
    weight = np.zeros_like(o)
    weight[valid] = io[valid] * ir[valid] / (io[valid] + ir[valid])
    residual = np.where(valid, o - r, 0.0)
    template = p[::-1, ::-1]
    rho = fftconvolve(weight * residual, template, mode="same")
    kappa = np.maximum(fftconvolve(weight, template**2, mode="same"), 0.0)
    support = fftconvolve(valid.astype(float), template**2, mode="same") / np.sum(p**2)
    accepted = (support >= minimum_support - 1e-8) & (kappa > 0)
    amplitude = np.full(o.shape, np.nan)
    sigma = np.full(o.shape, np.nan)
    snr = np.full(o.shape, np.nan)
    amplitude[accepted] = rho[accepted] / kappa[accepted]
    sigma[accepted] = 1.0 / np.sqrt(kappa[accepted])
    snr[accepted] = rho[accepted] / np.sqrt(kappa[accepted])
    return {"amplitude_native": amplitude, "sigma_white_native": sigma,
            "snr_white": snr, "kappa": kappa, "rho": rho,
            "support": support, "eligible": accepted}


def peaks(snr: np.ndarray, threshold: float = 7.0, distance_pix: int = 5) -> list[dict]:
    """Provisional pixel peaks only; NOT vetted astrophysical detections."""
    if not (math.isfinite(threshold) and threshold > 0 and distance_pix >= 1):
        raise ValueError("Invalid threshold or radius")
    arr = np.asarray(snr, float)
    clean = np.where(np.isfinite(arr), arr, -np.inf)
    local_max = maximum_filter(clean, size=2 * distance_pix + 1, mode="constant", cval=-np.inf)
    ys, xs = np.where((clean >= threshold) & (clean == local_max))
    return [{"y": int(y), "x": int(x), "snr_white": float(clean[y, x])}
            for y, x in sorted(zip(ys, xs), key=lambda yx: -clean[yx])]


def negative_tail_count(snr: np.ndarray, threshold: float = 7.0) -> int:
    """A diagnostic; -residual is NOT an independent jackknife."""
    a = np.asarray(snr, float)
    return int(np.count_nonzero(np.isfinite(a) & (a <= -threshold)))


def poisson_zero_rate_upper_95(area_time_deg2_days: float, *, validated: bool = False) -> float:
    """95% Poisson upper bound for zero detections given VALIDATED effective exposure.

    Exposure must already integrate injection completeness over map area/time/flux.
    Correlated observing visits, cuts and event duration belong in that integral.
    """
    if not validated:
        raise ValueError("Unvalidated footprint/completeness: no scientific rate limit")
    if not (math.isfinite(area_time_deg2_days) and area_time_deg2_days > 0):
        raise ValueError("Effective area-time must be finite and positive")
    return -math.log(0.05) / area_time_deg2_days
