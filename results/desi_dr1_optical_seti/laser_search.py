"""Search one DESI target and arm. Wavelengths are vacuum Angstroms in the SSB frame.

Flux density uses DESI units of 1e-17 erg/s/cm2/Angstrom. Integrated flux
uses 1e-17 erg/s/cm2. Positive input velocities recede in the relativistic
Doppler convention. Supply the complete target EXPID inventory separately.
"""
from dataclasses import dataclass, replace

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.sparse import dia_matrix
from scipy.stats import chi2

C_KMS = 299792.458
PARSEC_M = 3.085677581491367e16


@dataclass(frozen=True)
class Exposure:
    expid: int
    wave: np.ndarray
    flux: np.ndarray
    ivar: np.ndarray
    mask: np.ndarray
    resolution: np.ndarray
    template: np.ndarray  # stellar template convolved to this fibre/arm/exposure
    velocity_kms: float
    velocity_error_kms: float = 0.0  # independent RV error only

    @property
    def doppler(self):
        beta = self.velocity_kms / C_KMS
        return np.sqrt((1 + beta) / (1 - beta))


@dataclass(frozen=True)
class Hit:
    expid: int
    observed_A: float  # barycentric
    rest_A: float
    rest_error_A: float
    rest_fwhm_A: float  # Gaussian-equivalent, positive LSF moments
    flux17: float
    flux_error17: float
    snr: float
    shape_p: float
    spike_delta_chi2: float  # positive favours the spike alternative
    broad_delta_chi2: float  # positive favours the broadened alternative
    shape_pass: bool


@dataclass(frozen=True)
class Candidate:
    rest_A: float
    hits: tuple


@dataclass(frozen=True)
class Search:
    coverage: tuple  # intersection of eligible wavelength intervals, stellar rest A
    exposure_hits: dict  # all SNR triggers including failed shape checks
    candidates: tuple  # persistence and shape only; requires astrophysical/instrument vetting


@dataclass(frozen=True)
class Settings:
    snr: float = 7.
    shape_p_min: float = 0.001
    alternative_delta_chi2: float = 9.
    # Calibration floor per exposure, in native wavelength samples.
    calibration_pixels: float = 0.12
    centroid_sigmas: float = 3.


def convolve_template(resolution, intrinsic):
    """Intrinsic flux density must already be sampled on the native wave grid."""
    n = len(intrinsic)
    offsets = len(resolution) // 2 - np.arange(len(resolution))
    matrix = dia_matrix((np.asarray(resolution, dtype=float), offsets), shape=(n, n))
    return matrix @ np.asarray(intrinsic, dtype=float)


def _validate(e):
    w = np.asarray(e.wave)
    if w.ndim != 1 or not np.all(np.isfinite(w)) or not np.all(np.diff(w) > 0):
        raise ValueError("wave must be a finite strictly increasing grid")
    if any(np.shape(a) != w.shape for a in (e.flux, e.ivar, e.mask, e.template)):
        raise ValueError("flux, template, ivar and mask must match wave")
    r = np.asarray(e.resolution)
    if r.ndim != 2 or r.shape[1] != w.size or r.shape[0] < 3 or not r.shape[0] % 2:
        raise ValueError("resolution must be [odd ndiag, nwave]")
    if not np.isfinite(e.velocity_kms) or abs(e.velocity_kms) >= C_KMS:
        raise ValueError("invalid barycentric stellar velocity")
    if not np.isfinite(e.velocity_error_kms) or e.velocity_error_kms < 0:
        raise ValueError("invalid velocity error")


def _intervals(wave, centers):
    if not len(centers):
        return ()
    splits = np.split(centers, np.flatnonzero(np.diff(centers) != 1) + 1)
    return tuple((float(wave[g[0]]), float(wave[g[-1] + 1])) for g in splits)


def intersect_coverage(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a, b = left[i]
        c, d = right[j]
        if max(a, c) < min(b, d):
            out.append((max(a, c), min(b, d)))
        if b <= d:
            i += 1
        if d <= b:
            j += 1
    return tuple(out)


def _project(values, weights, x):
    """Remove a weighted affine continuum per window. x is weighted-centred."""
    mean = (values * weights).sum(1) / weights.sum(1)
    slope = (values * weights * x).sum(1) / (weights * x**2).sum(1)
    return values - mean[:, None] - slope[:, None] * x


def _alternative_gain(profile, y, weight, x, reference):
    p = _project(profile[None, :], weight[None, :], x[None, :])[0]
    norm = np.dot(weight, p**2)
    amplitude_score = max(np.dot(weight * y, p), 0.)
    return amplitude_score**2 / norm - reference if norm > 0 else -np.inf


def scan_exposure(e, settings=Settings()):
    """Profile the continuous adjacent-column mixture, affine baseline and amplitude.

    Only complete unmasked fitting windows are eligible. A 4-pixel sideband
    on each side constrains the affine residual in addition to template subtraction.
    Return threshold triggers and the complete list of eligible rest intervals.
    """
    _validate(e)
    w = np.asarray(e.wave, dtype=float)
    r = np.asarray(e.resolution, dtype=float)
    m, n = r.shape
    half = m // 2
    side = 4
    pad = half + side
    width = 2 * pad + 2  # includes support of the adjacent (j+1) column
    if n < width:
        return (), ()
    residual = np.asarray(e.flux, dtype=float) - np.asarray(e.template, dtype=float)
    ivar = np.asarray(e.ivar, dtype=float)
    good = (np.asarray(e.mask) == 0) & (ivar > 0) & np.isfinite(ivar) & np.isfinite(residual)
    valid = sliding_window_view(good, width).all(1)
    centers = np.flatnonzero(valid) + pad
    wave_width = np.gradient(w)
    # DESI DIA data[k,j] lands at output row j - (m//2-k).
    offsets = np.arange(-half, half + 1)
    good_center = np.ones(centers.size, dtype=bool)
    norms = []
    for j in (centers, centers + 1):
        col = r[:, j].T
        norm = (col * wave_width[j[:, None] + offsets]).sum(1)
        good_center &= np.isfinite(col).all(1) & (norm > 0) & np.isfinite(norm)
        norms.append(norm)
    j = centers[good_center]
    if not j.size:
        return (), ()
    coverage = _intervals(w / e.doppler, j)
    starts = j - pad
    weights = sliding_window_view(ivar, width)[starts]
    x = sliding_window_view(w, width)[starts] - w[j, None]
    x -= ((weights * x).sum(1) / weights.sum(1))[:, None]
    y = _project(sliding_window_view(residual, width)[starts], weights, x)
    p0 = np.zeros_like(y)
    p1 = np.zeros_like(y)
    p0[:, side:side+m] = r[:, j].T / norms[0][good_center, None]
    p1[:, side+1:side+m+1] = r[:, j+1].T / norms[1][good_center, None]
    u = _project(p0, weights, x)
    v = _project(p1 - p0, weights, x)
    a = (weights * y * u).sum(1)
    b = (weights * y * v).sum(1)
    c = (weights * u**2).sum(1)
    d = (weights * u * v).sum(1)
    f = (weights * v**2).sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        stationary = np.clip((a*d - b*c) / (b*d - a*f), 0., 1.)
        phases = np.column_stack((np.zeros(j.size), np.nan_to_num(stationary), np.ones(j.size)))
        norms2 = c[:, None] + 2*phases*d[:, None] + phases**2*f[:, None]
        scores = (a[:, None] + b[:, None]*phases) / np.sqrt(norms2)
    scores[~np.isfinite(scores)] = -np.inf
    choice = scores.argmax(1)
    t = phases[np.arange(j.size), choice]
    snr = scores[np.arange(j.size), choice]
    # Neighbouring intervals often identify the very same peak at an endpoint.
    # Suppress overlapping solutions in descending significance within one pixel.
    selected = []
    used = set()
    triggers = np.flatnonzero(snr >= settings.snr)
    for k in triggers[np.argsort(-snr[triggers], kind="stable")]:
        nearest = int(np.floor(j[k] + t[k] + 0.5))
        if nearest in used:
            continue
        # Suppress a monotonic shoulder even when it extends several pixels.
        left = k > 0 and j[k-1] == j[k]-1 and snr[k-1] > snr[k]
        right = k+1 < j.size and j[k+1] == j[k]+1 and snr[k+1] > snr[k]
        if left or right:
            continue
        used.update((nearest-1, nearest, nearest+1))
        selected.append(k)
    hits = []
    for k in sorted(selected):
        phase = t[k]
        p = u[k] + phase*v[k]
        profile = p0[k] + phase*(p1[k] - p0[k])
        norm = np.dot(weights[k], p**2)
        error = 1 / np.sqrt(norm)
        amplitude = snr[k] * error
        chi = np.dot(weights[k], (y[k] - amplitude*p)**2)
        p_shape = chi2.sf(chi, width - 4)
        pix = w[j[k]+1] - w[j[k]]
        observed = w[j[k]] + phase * pix
        w_window = w[starts[k]:starts[k]+width]
        mass = np.maximum(profile, 0.) * wave_width[starts[k]:starts[k]+width]
        mean = np.average(w_window, weights=mass)
        sigma = np.sqrt(np.average((w_window - mean)**2, weights=mass))
        derivative_info = f[k] - np.dot(weights[k] * v[k], p)**2 / norm
        if not np.isfinite(derivative_info) or derivative_info <= 0:
            continue
        stat_err = pix / (amplitude * np.sqrt(derivative_info))
        rest = observed / e.doppler
        # Exact derivative of ln(Doppler) with respect to relativistic beta.
        rv_err = rest * e.velocity_error_kms / (C_KMS * (1 - (e.velocity_kms/C_KMS)**2))
        rest_err = np.hypot(np.hypot(stat_err, settings.calibration_pixels*pix) / e.doppler, rv_err)
        reference = snr[k]**2
        # Compare positive single-pixel glitches at both candidate centre pixels.
        spike_gain = -np.inf
        for i in (pad, pad+1):
            spike = np.zeros(width)
            spike[i] = 1.
            spike_gain = max(spike_gain, _alternative_gain(spike, y[k], weights[k], x[k], reference))
        # Broaden by another instrumental sigma, giving ~sqrt(2) times the FWHM.
        ii = np.arange(width) - width // 2
        kernel = np.exp(-0.5 * (ii / max(sigma/pix, 0.1))**2)
        full = np.convolve(profile, kernel / kernel.sum())
        broad = full[width//2:width//2 + width]
        broad_gain = _alternative_gain(broad, y[k], weights[k], x[k], reference)
        shape_pass = p_shape >= settings.shape_p_min and max(spike_gain, broad_gain) <= settings.alternative_delta_chi2
        hits.append(Hit(int(e.expid), observed, rest, rest_err, 2.35482*sigma/e.doppler,
                        amplitude, error, snr[k], p_shape, spike_gain, broad_gain, bool(shape_pass)))
    return tuple(hits), coverage


def _compatible(h, g, settings):
    max_separation = min(settings.centroid_sigmas * np.hypot(h.rest_error_A, g.rest_error_A),
                         0.5 * min(h.rest_fwhm_A, g.rest_fwhm_A))
    return abs(h.rest_A - g.rest_A) <= max_separation


def scan_target(exposures, expected_expids, settings=Settings()):
    """Search one camera arm across every distinct exposure in the input inventory."""
    exposures = tuple(exposures)
    expected = tuple(expected_expids)
    ids = [int(e.expid) for e in exposures]
    if len(expected) < 2 or len(set(expected)) != len(expected):
        raise ValueError("expected_expids must contain at least two distinct EXPIDs")
    if len(ids) != len(set(ids)) or set(ids) != set(expected):
        raise ValueError("incomplete or duplicate exposure records for this target and arm")
    all_hits = {}
    coverage = None
    for e in exposures:
        hits, eligible = scan_exposure(e, settings)
        all_hits[e.expid] = hits
        coverage = eligible if coverage is None else intersect_coverage(coverage, eligible)
    passing = [[h for h in hits if h.shape_pass] for hits in all_hits.values()]
    passing.sort(key=len)
    candidates = []
    for anchor in passing[0]:
        matched = [anchor]
        for others in passing[1:]:
            nearby = [h for h in others if all(_compatible(h, g, settings) for g in matched)]
            if not nearby:
                break
            matched.append(min(nearby, key=lambda h: abs(h.rest_A-anchor.rest_A)))
        if len(matched) == len(exposures):
            rest = np.average([h.rest_A for h in matched], weights=[1/h.rest_error_A**2 for h in matched])
            if any(lo <= rest <= hi for lo, hi in coverage):
                candidates.append(Candidate(float(rest), tuple(sorted(matched, key=lambda h: h.expid))))
    return Search(coverage, all_hits, tuple(candidates))


def inject_line(e, rest_A, flux17):
    """Add a flux-conserving unresolved signal at a continuous subpixel phase."""
    if not np.isfinite(flux17) or flux17 <= 0:
        raise ValueError("injected integrated flux must be positive and finite")
    w = np.asarray(e.wave, dtype=float)
    observed = rest_A * e.doppler
    j = np.searchsorted(w, observed, side="right") - 1
    half = e.resolution.shape[0] // 2
    if j-half < 0 or j+1+half >= w.size:
        raise ValueError("injection outside arm interior")
    phase = (observed - w[j]) / (w[j+1] - w[j])
    flux = np.asarray(e.flux, dtype=float).copy()
    widths = np.gradient(w)
    for col, fraction in ((j, 1-phase), (j+1, phase)):
        rows = slice(col-half, col+half+1)
        response = np.asarray(e.resolution[:, col], dtype=float)
        integral = response @ widths[rows]
        if not np.isfinite(integral) or integral <= 0:
            raise ValueError("invalid resolution column")
        flux[rows] += flux17 * fraction * response / integral
    return replace(e, flux=flux)


def isotropic_power_w(flux17, distance_pc):
    """Isotropic-equivalent, exposure-averaged line power in W, no extinction correction."""
    if not all(np.isfinite(v) and v > 0 for v in (flux17, distance_pc)):
        raise ValueError("flux and distance must be positive and finite")
    return 4 * np.pi * (PARSEC_M*distance_pc)**2 * flux17 * 1e-20
