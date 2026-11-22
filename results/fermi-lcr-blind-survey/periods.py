"""Blind LCR periodicity scan and catalog-wide Monte Carlo calibration.

Inputs are fixed-index photon flux detections, upper limits, and MJD times.
Surrogates describe the measured flux distribution conditional on the window.
The continuum is a fitted model, so the reported FAPs are model-conditional.
"""
import hashlib
import math

import numpy as np
from scipy.stats import beta

PILOT = 96
REPLICATES = 1023
INJECTIONS = 16
AMPLITUDES = (0.25, 0.50, 0.90)
INJECTION_PERIODS = (240., 540., 1080.)
ALPHA = 0.05


def _column(curve, key, fallback):
    if key in curve:
        return np.asarray(curve[key], dtype=float)
    return np.asarray(fallback, dtype=float)


def _points(curve):
    t = _column(curve, "t", curve.get("time", []))
    flux = _column(curve, "flux", np.full_like(t, np.nan))
    ul = _column(curve, "ul", curve.get("upper", np.full_like(t, np.nan)))
    ok = np.isfinite(t)
    flux = np.where(np.isfinite(flux) & (flux > 0), flux, np.nan)
    ul = np.where(np.isfinite(ul) & (ul > 0), ul, np.nan)
    order = np.argsort(t[ok])
    return t[ok][order], flux[ok][order], ul[ok][order]


def _grid(t, cadence):
    """Infer missing integer bins from the native archive sampling."""
    gaps = np.diff(t)
    low, high, max_bins = (26., 33., 350) if cadence == "monthly" else (6.5, 7.5, 1500)
    adjacent = gaps[(gaps >= low) & (gaps <= high)]
    if not adjacent.size:
        return None
    step = np.median(adjacent)
    indices = np.r_[0, np.cumsum(np.rint(gaps / step).astype(int))]
    if np.any(np.diff(indices) < 1) or indices[-1] > max_bins:
        return None
    step, origin = np.polyfit(indices, t, 1)
    if np.max(np.abs(t - origin - indices * step)) > 0.20 * step:
        return None
    return indices, float(step), float(origin)


def _continuum(grid, mask):
    """Profile the window-convolved expected periodogram over spectral shapes.

    The kernel accounts for the detection mask and subtraction of its sample
    mean. The underlying stationary process is four observation spans long.
    """
    n = len(grid)
    length = 4 * n
    window = np.fft.fft(mask.astype(float), n=length)
    harmonics = np.arange(4, (n // 2) * 4 + 1, 4)
    j = np.arange(length)
    window_kernel = np.abs(window[(harmonics[:, None] - j) % length]
                           - window[harmonics, None] * window[-j] / mask.sum()) ** 2
    observed = np.abs(np.fft.rfft(grid, n=length)[harmonics]) ** 2
    freq = np.maximum(np.abs(np.fft.fftfreq(length)), 1. / length)
    slopes = np.arange(0., 3.26, 0.25)
    # White power relative to the red component at Nyquist.
    floors = np.array([0., 0.3, 1., 3., 10., 30., 100., 300., 1000.])
    red = (freq[:, None] / 0.5) ** (-slopes)
    red[0] = 0.
    convolved = (window_kernel @ red)[:, :, None] + window_kernel[:, 1:].sum(axis=1)[:, None, None] * floors
    convolved = np.maximum(convolved, 1.e-14)
    scale = (observed[:, None, None] / convolved).mean(axis=0)
    cost = np.log(convolved).mean(axis=0) + np.log(np.maximum(scale, 1.e-30))
    i, j = np.unravel_index(np.argmin(cost), cost.shape)
    model = red[:length // 2 + 1, i] + floors[j]
    model[0] = 0.
    return model / np.max(model), float(slopes[i]), float(floors[j])


def _remap(values, ordered):
    mapped = np.empty_like(values)
    indices = np.argsort(values, axis=1)
    np.put_along_axis(mapped, indices, np.broadcast_to(ordered, values.shape), axis=1)
    return mapped


class Surrogates:
    def __init__(self, data, mask, spectrum, rng):
        self.mask = mask
        self.n = mask.size
        self.length = 4 * self.n
        self.observed = np.sort(data)
        self.marginal = np.quantile(data, (np.arange(self.length) + 0.5) / self.length)
        self.amplitudes = np.sqrt(spectrum)
        self.rng = rng

    def sample(self, count):
        # Random Fourier amplitudes retain intrinsic periodogram scatter. Using
        # fixed observed Fourier magnitudes would preserve a real signal peak.
        shape = (count, len(self.amplitudes))
        z = (self.rng.normal(size=shape) + 1j * self.rng.normal(size=shape)) * self.amplitudes
        z[:, -1] = self.rng.normal(size=count) * self.amplitudes[-1] * np.sqrt(2.)
        target = np.abs(z)
        values = np.fft.irfft(z, n=self.length, axis=1)
        for _ in range(8):
            values = _remap(values, self.marginal)
            z = np.fft.rfft(values, axis=1)
            z *= target / np.maximum(np.abs(z), 1.e-40)
            values = np.fft.irfft(z, n=self.length, axis=1)
        # Same censoring window; enforce EXACT measured marginal in that
        # window. Retained-spectrum diagnostics below include this last step.
        return _remap(values[:, :self.n][:, self.mask], self.observed)


def _fourier(values, mask, split, freq_bins, step):
    """Oversampled, mean-subtracted, unweighted Fourier coefficients."""
    values = np.atleast_2d(values)
    n = mask.size
    indices = np.flatnonzero(mask)
    if split == 0:
        use = np.ones(len(indices), dtype=bool)
    elif split == 1:
        use = indices < n / 2
    else:
        use = indices >= n / 2
    selected = values[:, use]
    padded = np.zeros((len(values), n))
    padded[:, indices[use]] = selected - selected.mean(axis=1, keepdims=True)
    return np.fft.rfft(padded, n=4*n, axis=1)[:, freq_bins]


def _scan(values, mask, bins, norms, step):
    coeff = [_fourier(values, mask, split, bins, step) for split in range(3)]
    powers = [np.abs(c) ** 2 / norm for c, norm in zip(coeff, norms)]
    best = [p.argmax(axis=1) for p in powers]
    rows = np.arange(len(best[0]))
    peak = best[0]
    # Each disjoint half must peak within one half-span Fourier resolution.
    proximity = [np.abs(bins[x] - bins[peak]) <= 8 for x in best[1:]]
    phase_gap = np.abs(np.angle(coeff[1][rows, peak] * coeff[2][rows, peak].conj()))
    stable = proximity[0] & proximity[1] & (phase_gap <= np.pi / 3.)
    return (powers[0][rows, peak], bins[peak], stable, phase_gap,
            bins[best[1]], bins[best[2]])


def _coefficients(t, freq):
    phase = 2 * np.pi * freq * t
    return np.column_stack((np.ones(len(t)), np.cos(phase), np.sin(phase)))


def _direction(t, flux, freq):
    ok = np.isfinite(flux)
    if ok.sum() < 12:
        return None
    x = _coefficients(t[ok], freq)
    if np.linalg.cond(x) > 10:
        return None
    fit = np.linalg.lstsq(x, flux[ok], rcond=None)[0][1:]
    return fit / max(np.linalg.norm(fit), 1.e-50)


def _censor_margin(t, flux, ul, freq, direction):
    # Lower bound on the signed coherent amplitude in the fitted direction,
    # allowing each published upper limit to vary separately from zero to UL.
    keep = np.isfinite(flux) | np.isfinite(ul)
    x = _coefficients(t[keep], freq)
    if len(x) < 12 or np.linalg.cond(x) > 10:
        return None
    weights = direction @ np.linalg.pinv(x)[1:]
    flux, ul = flux[keep], ul[keep]
    detected = np.isfinite(flux)
    bound = weights[detected] @ flux[detected] + np.minimum(weights[~detected], 0.) @ ul[~detected]
    return float(bound)


def _is_alias(frequency, span):
    return np.any(np.abs(np.asarray(frequency)[..., None] - np.array([1., 2.]) / 365.25) <= 1. / span, axis=-1)


def _coverage(t, flux, midpoint, span):
    observed = t[np.isfinite(flux)]
    return bool(all(((observed >= start) & (observed < start + span/4)).sum() >= 3
                    for start in midpoint + span * np.array([-.5, -.25, 0., .25])))


def analyze(source, cadences):
    """Return JSON-safe source row, float32 null maxima, and injection batches."""
    # Select solely by native product availability, before testing coverage or power.
    primary = next((c for c in ("monthly", "weekly") if c in cadences), None)
    row = dict(source=str(source), primary_cadence=primary, status="monthly and weekly unavailable",
               passes=False, period_days=None, p_global=None)
    empty = (row, np.array([], dtype=np.float32), [])
    if primary is None:
        return empty
    t, flux, ul = _points(cadences[primary])
    row[primary + "_input"] = len(t)
    row[primary + "_detections"] = int(np.isfinite(flux).sum())
    if len(t) < 48:
        row['status'] = f"fewer than 48 {primary} timestamps"
        return empty
    axis = _grid(t, primary)
    if axis is None:
        row['status'] = f"irregular or duplicate {primary} grid"
        return empty
    idx, step, origin = axis
    n = int(idx[-1]) + 1
    mask = np.zeros(n, dtype=bool)
    mask[idx] = np.isfinite(flux)
    values_raw = flux[np.isfinite(flux)]
    span = n * step
    midpoint = origin + (n - 1) * step / 2.
    row.update(start_mjd=float(t[0]), end_mjd=float(t[-1]), detection_fraction=float(mask.mean()))
    row[primary + "_step_days"] = step
    row[primary + "_grid"] = n
    if len(values_raw) < 36 or min(mask[:(n+1)//2].mean(), mask[(n+1)//2:].mean()) < 0.20 or not _coverage(t, flux, midpoint, span):
        row['status'] = "sparse detections or quarter coverage"
        return empty
    rms = float(np.std(values_raw))
    mean = float(np.mean(values_raw))
    if not np.isfinite(rms) or rms <= 0.:
        row['status'] = "zero or invalid flux variance"
        return empty
    rng = np.random.default_rng(int.from_bytes(hashlib.sha256(("lcr-blind-v1|" + str(source)).encode()).digest()[:8], "little"))
    values = (values_raw - mean) / rms
    grid = np.zeros(n)
    grid[mask] = values
    spectrum, slope, white = _continuum(grid, mask)
    sampler = Surrogates(values, mask, spectrum, rng)
    # Frequencies are fixed before looking at powers. Four observed cycles
    # minimum and a roughly half-year lower-period limit in both products.
    minimum_bins = 6 if primary == "monthly" else 26
    bins = np.arange(16, int(4*n/minimum_bins) + 1)
    if len(bins) < 8:
        row['status'] = "period range too short"
        return empty
    pilot = sampler.sample(PILOT)
    norms = []
    for split in range(3):
        powers = np.abs(_fourier(pilot, mask, split, bins, step)) ** 2
        norms.append(np.maximum(powers.mean(axis=0), 1.e-10))
    observed = _scan(values, mask, bins, norms, step)
    score, k, stable, phase, left, right = [x[0] for x in observed]
    period = 4*span / k
    frequency = 1./period
    # Independent-of-scoring diagnostics compare broad spectral bands across
    # the entire native frequency range. Use medians to resist a narrow line.
    check_bins = np.arange(4, n//2*4+1, 4)
    observed_shape = np.abs(_fourier(values, mask, 0, check_bins, step)[0]) ** 2
    pilot_shape = np.abs(_fourier(pilot, mask, 0, check_bins, step)) ** 2
    bands = np.array_split(np.arange(len(check_bins)), 5)
    band_ratios = []
    model_ok = True
    for band in bands:
        sims = np.median(pilot_shape[:, band], axis=1)
        obs = np.median(observed_shape[band])
        band_ratios.append(float(obs / max(np.median(sims), 1.e-15)))
        # A known measurement window is reapplied to every realization. Flag
        # extreme spectral mismatch instead of trusting an inadequate null.
        model_ok &= sims.min() <= obs <= sims.max()
    direction = _direction(t-origin, flux, frequency)
    bounds = []
    confirms = {}
    if direction is not None:
        for side in (t < midpoint, t >= midpoint):
            bounds.append(_censor_margin(t[side]-origin, flux[side], ul[side], frequency, direction))
        confirmations = ("weekly", "3day") if primary == "monthly" else ("3day",)
        for cadence in confirmations:
            if cadence not in cadences:
                continue
            ct, cf, cu = _points(cadences[cadence])
            if not _coverage(ct, cf, midpoint, span):
                continue
            orient = _direction(ct-origin, cf, frequency)
            if orient is None:
                continue
            offset = float(np.degrees(np.arccos(np.clip(direction @ orient, -1., 1.))))
            confirms[cadence] = offset
    null = []
    for count in range(0, REPLICATES, 32):
        surrogate = sampler.sample(min(32, REPLICATES-count))
        powers = np.abs(_fourier(surrogate, mask, 0, bins, step)) ** 2 / norms[0]
        null.extend(powers.max(axis=1))
    row.update(status="tested", period_days=float(period), score=float(score),
               p_source_scan=(int(np.sum(np.array(null) >= score))+1)/(REPLICATES+1),
               pmin_days=float(4*span/bins[-1]), pmax_days=float(4*span/bins[0]),
               n_frequencies=len(bins), psd_slope=slope, white_to_red_nyquist=white,
               spectral_band_ratios=band_ratios, spectrum_ok=bool(model_ok),
               left_period_days=float(4*span/left), right_period_days=float(4*span/right),
               half_phase_gap_deg=float(np.degrees(phase)), halves_ok=bool(stable),
               censor_half_amplitude_lower=bounds, censor_ok=len(bounds)==2 and all(b is not None and b>0 for b in bounds),
               annual_alias=bool(_is_alias(frequency, span)), cadence_phase_gaps_deg=confirms,
               cadence_ok=bool(confirms) and all(x <= 60. for x in confirms.values()))
    injections = []
    for ip in INJECTION_PERIODS:
        if not row['pmin_days'] <= ip <= row['pmax_days']:
            continue
        for a in AMPLITUDES:
            base = sampler.sample(INJECTIONS) * rms + mean
            phase = 2*np.pi*rng.random((INJECTIONS, 1))
            signal = 1. + a * np.cos(2*np.pi*step*np.flatnonzero(mask)/ip + phase)
            simulated = (base * signal - mean) / rms
            # Apply the same variance convention as the reference distribution.
            simulated /= np.maximum(simulated.std(axis=1, keepdims=True), 1.e-12)
            simulated_scores, peaks, stable, _, _, _ = _scan(simulated, mask, bins, norms, step)
            detected_freq = peaks/(4*span)
            good = (np.abs(detected_freq - 1./ip) <= 1./span) & stable & ~_is_alias(detected_freq, span)
            injections.append((ip, a, np.where(good, simulated_scores, -1.).astype(np.float32)))
    return row, np.asarray(null, dtype=np.float64), injections


def finalize(records):
    """Calibrate source scores using catalog maxima; recover injected signals."""
    rows = [item[0] for item in records]
    eligible = [r for r in records if len(r[1])]
    injection_rows = []
    summary = dict(roster=len(rows), tested=len(eligible), replicates=REPLICATES, alpha=ALPHA,
                   pilot=PILOT, accepted=0, global_threshold=None, aggregate_recovery=[])
    if not eligible:
        return rows, injection_rows, summary
    maxima = np.max(np.stack([r[1] for r in eligible]), axis=0)
    # For reporting injections use the exact same conservative Monte Carlo
    # gate (95% one-sided binomial upper bound AND corrected rank p <= .05).
    limit = max((k for k in range(REPLICATES)
                 if beta.ppf(.95, k+1, REPLICATES-k) <= ALPHA
                 and (k+1)/(REPLICATES+1) <= ALPHA), default=-1)
    threshold = float(np.sort(maxima)[-(limit+1)]) if limit >= 0 else math.inf
    summary['global_threshold'] = threshold
    # Retain the calibration distribution for independent verification.
    summary['catalog_null_maxima'] = maxima.astype(float).tolist()
    summary['global_exceedance_limit'] = limit
    aggregate = {}
    for row, _, injected in eligible:
        exceed = int(np.sum(maxima >= row['score']))
        row['p_global'] = (1+exceed)/(1+REPLICATES)
        row['global_exceedances'] = exceed
        row['p_global_mc_upper95'] = float(beta.ppf(.95, exceed+1, REPLICATES-exceed)) if exceed < REPLICATES else 1.
        row['global_ok'] = exceed <= limit
        row['passes'] = bool(row['global_ok'] and row['halves_ok'] and row['cadence_ok'] and row['spectrum_ok']
                             and row['censor_ok'] and not row['annual_alias'])
        summary['accepted'] += row['passes']
        for ip, a, scores in injected:
            recovered = int(np.sum(np.asarray(scores) > threshold))
            injection_rows.append(dict(source=row['source'], primary_cadence=row['primary_cadence'],
                                       period_days=ip, semi_amplitude=a,
                                       n=len(scores), recovered=recovered))
            totals = aggregate.setdefault((ip, a), [0, 0])
            totals[0] += recovered
            totals[1] += len(scores)
    for (ip, a), (k, n) in sorted(aggregate.items()):
        summary['aggregate_recovery'].append(dict(period_days=ip, semi_amplitude=a,
                                        recovered=k, n=n, fraction=k/n))
    return rows, injection_rows, summary


METHOD = """The blind scan used native monthly fixed-index photon-flux detections whenever that product was available. Sources without an acquired monthly file used their native weekly product if available. This choice preceded all quality and power tests. Each source entered at most once. Each source needed at least 48 timestamps, 36 detections, 20% detection coverage in each half and three detections per time quarter. The native grid had to fit within 0.20 bin. The unweighted, mean-subtracted Fourier scan used four samples per independent Fourier interval and periods between six monthly bins (26 weekly bins for the fallback) and one-quarter of the baseline.

We fit a red power law plus a white floor by profiled window-convolved Whittle likelihood over native harmonics. The stochastic process extended four spans, with its lowest nonzero frequency at one-quarter of the inverse observed span. Each realization used independent Gaussian Fourier coefficients, eight iterative spectral-amplitude/rank adjustments, a four-to-one crop, the exact detection mask and exact detected-flux ranks. A separate 96-surrogate pilot calibrated frequency-dependent score normalization. We retained ratios of observed to simulated median powers in five frequency bands and flagged any band outside the pilot range. These are plug-in inferences conditional on the fitted stationary spectrum, detected-flux distribution and sampling mask. The procedure covers measured flux variability and does not reconstruct the latent flux below detection thresholds.

Each of 1,023 synthetic catalogs takes the maximum normalized power over every period and every tested source. The catalog simulation assumes independent source noise; shared exposure systematics and fitted continuum uncertainty are outside that calibration. The reported family-wise p is (exceedances+1)/1,024; discovery screening also requires its one-sided 95% binomial upper limit to be at most 0.05. Both nonoverlapping halves must independently peak within one half-span Fourier interval of the full-series peak and agree in phase within 60 degrees. The fitted harmonic in each half must retain positive projected amplitude for every assignment of individual upper-limit bins between zero and their reported bound. These individual bounds are sensitivity scenarios, not simultaneous confidence intervals. All usable short-cadence products must agree in fixed-period phase within 60 degrees, with at least one available. For weekly primaries, only the 3-day product can confirm the phase. These products share photons with the primary series. Frequencies within one full-span Fourier interval of the annual or semiannual cycle are flagged.

Injection recovery used 16 independent source-specific realizations per supported source, period and amplitude with random phase. A sinusoid multiplied positive simulated fluxes at the observed detection epochs by 1+a*cos(phase). The recovered peak must clear the same survey threshold, match frequency within one full-span bin, pass the disjoint-half test and avoid the annual veto. Recovery measures the conditional primary-series timing sensitivity. It holds the measured window and fitted continuum fixed and does not include censoring changes, short-cadence confirmation, counterpart identification or literature screening. Descriptive pooled fractions weight eligible sources equally at each injected period and amplitude."""
