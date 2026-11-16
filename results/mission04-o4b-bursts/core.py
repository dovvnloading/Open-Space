"""Two-site burst search. Times are measured from a common integer GPS origin.

The statistic assumes waveforms approximately related by a delay and a real
scale factor. The scan covers either sign. Noise calibration is empirical.
"""
from dataclasses import dataclass

import numpy as np
from scipy import fft, ndimage, signal

RATE = 2048
DURATIONS = (0.125, 0.5, 1.0)
BANDS = ((32, 96), (96, 224), (224, 448), (448, 896))
LAGS = np.arange(-20, 21) * 0.0005
EXCESS_GATE = 1.5
HALF = 16
MIN_TRAIN = 8
HIGHPASS = signal.firwin(513, 24, pass_zero=False, fs=RATE)


@dataclass
class Scale:
    duration: float
    freqs: np.ndarray
    bands: list
    phases: list
    spectra: np.ndarray
    excess: np.ndarray
    reference: dict


@dataclass
class Prepared:
    first_second: int
    live: np.ndarray
    scales: list


@dataclass(frozen=True)
class Trigger:
    time: float
    score: float
    lag_seconds: float
    band: tuple
    duration: float
    coherence: float
    sign: int


@dataclass
class Slide:
    live_seconds: int
    triggers: list
    nominal_tiles: int
    gated_tiles: int


def _condition(strain, fs, good):
    raw = np.asarray(strain, dtype=np.float64)
    if fs not in (2048, 4096, 16384) or raw.ndim != 1 or raw.size != fs * len(good):
        raise ValueError("supply complete seconds at 2048, 4096, or 16384 Hz")
    finite = np.all(np.isfinite(raw.reshape(-1, fs)), axis=1)
    valid = finite & good
    # One whole adjacent second on both sides covers the longest half-window
    # and the combined compact support of the resampling and high-pass filters.
    live = ndimage.binary_erosion(valid, structure=np.ones(3, dtype=bool))
    raw = np.where(valid.repeat(fs), raw, 0.0)
    if fs > RATE:
        raw = signal.resample_poly(raw, 1, fs // RATE)
    cleaned = signal.oaconvolve(raw, HIGHPASS, mode="same")
    return cleaned, live


def _build_scale(cleaned, eligible, first, duration, ref):
    ns = len(eligible)
    size = round(RATE * duration)
    hop = size // 2
    frame_count = RATE // hop
    grid = fft.rfftfreq(size, 1 / RATE)
    # Mask the mains harmonics plus at least one resolution bin either side.
    line_distance = np.abs(grid - 60 * np.rint(grid / 60))
    keep = (grid >= BANDS[0][0]) & (grid < BANDS[-1][1]) & (line_distance > max(1, RATE / size))
    freqs = grid[keep]
    bands = [np.flatnonzero((freqs >= low) & (freqs < high)) for low, high in BANDS]
    phases = [np.exp(-2j * np.pi * np.outer(freqs[indices], LAGS)).astype(np.complex64)
              for indices in bands]
    spectra = np.zeros((ns, frame_count, len(freqs)), dtype=np.complex64)
    excess = np.full((ns, frame_count, len(BANDS)), -np.inf, dtype=np.float32)
    # Compute the small raw-FFT batches in double precision (physical strain
    # powers can be far below the normal float32 range).
    padded = np.pad(cleaned, (size // 2, size // 2))
    windows = np.lib.stride_tricks.sliding_window_view(padded, size)[hop // 2::hop]
    windows = windows[:ns * frame_count].reshape(ns, frame_count, size)
    taper = signal.windows.hann(size, sym=False)
    reference = {} if ref is None else ref.reference
    valid_psd = np.zeros(ns, dtype=bool)
    gps_seconds = first + np.arange(ns)
    block_ids = gps_seconds // (2 * HALF)
    for block in range(int(block_ids[0]), int(block_ids[-1]) + 1):
        left = max(0, 2 * HALF * block - first)
        right = min(ns, 2 * HALF * (block + 1) - first)
        if not eligible[left:right].any():
            continue
        raw = fft.rfft(windows[left:right] * taper, axis=-1)[..., keep]
        local_epoch = gps_seconds[left:right]
        for half_id in (2 * block, 2 * block + 1):
            target = (local_epoch // HALF == half_id) & eligible[left:right]
            if not target.any():
                continue
            if ref is None:
                # Opposite half of this 32-s epoch, trimming both its edges.
                # A target cannot raise its own short-burst noise estimate.
                train = ((local_epoch // HALF == (half_id ^ 1))
                         & (local_epoch % HALF > 0) & (local_epoch % HALF < HALF - 1)
                         & eligible[left:right])
                if train.sum() < MIN_TRAIN:
                    continue
                power = np.abs(raw[train]) ** 2
                estimate = np.median(power, axis=(0, 1)) / np.log(2)
                if not (np.all(np.isfinite(estimate)) and np.all(estimate > 0)):
                    continue
                reference[half_id] = np.sqrt(estimate)
            if half_id not in reference:
                continue
            normed = raw[target] / reference[half_id]
            if not np.all(np.isfinite(normed)) or np.abs(normed).max() > 1e17:
                continue
            selected = np.flatnonzero(target) + left
            valid_psd[selected] = True
            spectra[selected] = normed
            norm_power = np.abs(normed) ** 2
            for band_id, indices in enumerate(bands):
                energy = np.sum(norm_power[..., indices], axis=-1)
                excess[selected, :, band_id] = (energy - len(indices)) / np.sqrt(2 * len(indices))
    return Scale(duration, freqs, bands, phases, spectra, excess, reference), valid_psd


def prepare(strain, fs, good_seconds, *, reference=None, first_second=0):
    """Condition one detector and make tiles. good_seconds must already include all vetoes.

    first_second is relative to the common experiment origin (zero for a single
    full frame). For injected short crops pass the original detector Prepared
    as reference and preserve the original integer origin via first_second.
    All three tile resolutions must be valid to earn a second of exposure.
    """
    good = np.asarray(good_seconds)
    if good.ndim != 1 or not good.size or not np.all(np.isin(good, (False, True))):
        raise ValueError("good_seconds must be a non-empty 1-D binary array")
    if int(first_second) != first_second:
        raise ValueError("first_second must be integral")
    first = int(first_second)
    cleaned, eligible = _condition(strain, fs, good.astype(bool))
    scales = []
    live = eligible.copy()
    for i, duration in enumerate(DURATIONS):
        scale, ok = _build_scale(cleaned, eligible, first, duration,
                                 None if reference is None else reference.scales[i])
        scales.append(scale)
        live &= ok
    return Prepared(first, live, scales)


def _run_slide(h, l, shift):
    ns = h.live.size
    first = max(0, -shift)
    end = min(ns, ns - shift)
    allowed = np.zeros(ns, dtype=bool)
    if first < end:
        allowed[first:end] = h.live[first:end] & l.live[first + shift:end + shift]
    nlive = int(allowed.sum())
    best = {}
    gated = 0
    nominal = 0
    for hs, ls in zip(h.scales, l.scales):
        frames = hs.spectra.shape[1]
        nominal += nlive * frames * len(BANDS)
        if not nlive:
            continue
        for band_id, indices in enumerate(hs.bands):
            zh = hs.excess[first:end, :, band_id]
            zl = ls.excess[first + shift:end + shift, :, band_id]
            sec, frame = np.nonzero(allowed[first:end, None] & (zh >= EXCESS_GATE) & (zl >= EXCESS_GATE))
            gated += len(sec)
            # Limit temporary memory even during a prolonged disturbance.
            for start in range(0, len(sec), 256):
                pick = slice(start, start + 256)
                a, f = sec[pick], frame[pick]
                hspec = hs.spectra[first + a, f][:, indices]
                lspec = ls.spectra[first + a + shift, f][:, indices]
                # conj(L) holds +2*pi*f*delay for a delayed copy at L.
                cross = hspec * lspec.conj()
                phase = hs.phases[band_id]
                # No threaded BLAS launch for each small gated tile batch.
                corr = np.einsum("tf,fl->tl", cross, phase, optimize=False).real
                peak = np.argmax(np.abs(corr), axis=1)
                proj = corr[np.arange(len(peak)), peak]
                eh = zh[a, f].astype(np.float64) * np.sqrt(2 * len(indices)) + len(indices)
                el = zl[a, f].astype(np.float64) * np.sqrt(2 * len(indices)) + len(indices)
                coherence = np.clip(np.abs(proj) / np.sqrt(eh * el), 0, 1)
                balance = 2 * np.sqrt(eh * el) / (eh + el)
                scores = np.minimum(zh[a, f], zl[a, f]) * coherence * balance
                for i, score in enumerate(scores):
                    second = int(first + a[i])
                    if second in best and score <= best[second].score:
                        continue
                    best[second] = Trigger(
                        time=h.first_second + second + (int(f[i]) + 0.5) / frames,
                        score=float(score), lag_seconds=float(LAGS[peak[i]]), band=BANDS[band_id],
                        duration=hs.duration, coherence=float(coherence[i]), sign=1 if proj[i] >= 0 else -1,
                    )
    return Slide(nlive, [best[k] for k in sorted(best)], nominal, gated)


def scan(h, l, shifts=(0,)):
    """Return max over all searched tiles/lags/signs for each eligible second.

    A live second with no double excess-power gate crossing has score zero and
    is omitted from .triggers. +shift compares H(t) with L(t+shift); a positive
    residual lag means L arrives later. Slides do not wrap or share boundaries
    between separately prepared frames. Slides reuse detector data.
    """
    if h.first_second != l.first_second or len(h.live) != len(l.live):
        raise ValueError("both detectors must share the same time grid")
    shifts = tuple(shifts)
    if len(set(shifts)) != len(shifts) or any(int(s) != s or abs(s) == 1 for s in shifts):
        raise ValueError("use unique integer-second shifts with magnitude zero or >= 2")
    return {int(s): _run_slide(h, l, int(s)) for s in shifts}


@dataclass(frozen=True)
class Confirmation:
    peak_time: float
    lag_seconds: float
    sign: int
    correlation: float
    coherence: float
    signal_to_offsource: float


def confirm(h_strain, l_strain, fs, time, *, band=(32, 896), duration=0.25,
            first_second=0, good_h=None, good_l=None):
    """Independent raw-strain time-domain diagnostic around an already found trigger.

    The two series are independently bandpassed and whitened by low-order linear
    prediction trained off source. Direct summation of normalized products
    examines sliding windows, integer sample lags, and both polarities. The
    ratio to neighboring off-source products is a descriptive effect size,
    not a p-value. It must be calibrated separately if used for selection.
    Supply original veto masks when present. None indicates every second is
    originally usable. A 17-second *good* crop is generally sufficient.
    """
    arrays = [np.asarray(a, dtype=np.float64) for a in (h_strain, l_strain)]
    if any(a.ndim != 1 for a in arrays) or len(arrays[0]) != len(arrays[1]):
        raise ValueError("require aligned 1-D detector arrays")
    samples = len(arrays[0])
    if samples % fs or fs not in (2048, 4096, 16384) or not 0 < duration <= 1:
        raise ValueError("unsupported sampling grid or duration")
    if not (0 < band[0] < band[1] < fs / 2):
        raise ValueError("band must lie below Nyquist")
    second = int(np.floor(time - first_second))
    # Expand to 8 s either side, stopping before either detector's missing data.
    good = np.ones(samples // fs, dtype=bool)
    for a, mask in zip(arrays, (good_h, good_l)):
        if mask is not None:
            mask = np.asarray(mask)
            if mask.shape != good.shape or not np.all(np.isin(mask, (0, 1))):
                raise ValueError("invalid confirmation quality mask")
            good &= mask.astype(bool)
    left, right = max(0, second - 8), min(len(good), second + 9)
    if not left <= second < right:
        raise ValueError("candidate falls outside the input")
    for a in arrays:
        good[left:right] &= np.isfinite(a[left * fs:right * fs]).reshape(-1, fs).all(axis=1)
    if not good[second]:
        raise ValueError("candidate epoch fails data quality")
    before = np.flatnonzero(~good[left:second])
    after = np.flatnonzero(~good[second:right])
    if before.size:
        left += int(before[-1]) + 1
    if after.size:
        right = second + int(after[0])
    t0 = time - (first_second + left)
    n = (right - left) * fs
    local = np.arange(n) / fs
    train = (local > 0.5) & (local < (right - left - 0.5)) & (np.abs(local - t0) > 1.5)
    if np.count_nonzero(train) < 4 * fs or min(t0, right - left - t0) < 1.1:
        raise ValueError("insufficient good on/off-source conditioning data")
    # Reflection padding and half-second filter guards avoid using bad data.
    sos = signal.butter(4, band, fs=fs, btype="bandpass", output="sos")
    order = 12
    residuals = []
    from scipy.linalg import solve_toeplitz
    for a in arrays:
        segment = a[left * fs:right * fs]
        scale = np.median(np.abs(segment[train] - np.median(segment[train])))
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("degenerate off-source strain")
        filtered = signal.sosfiltfilt(sos, segment / scale)
        # Autocovariance pairs cannot cross disjoint chunks or enter on-source.
        cov = np.array([np.mean((filtered[k:] * filtered[:n-k])[train[k:] & train[:n-k]])
                        for k in range(order + 1)])
        if not np.isfinite(cov).all() or cov[0] <= 0:
            raise ValueError("degenerate bandpassed off-source strain")
        cov[0] *= 1.0001
        taps = solve_toeplitz(cov[:-1], cov[1:])
        residual = signal.lfilter(np.r_[1., -taps], [1.], filtered)
        scale = np.median(np.abs(residual[train])) / 0.67449
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("degenerate off-source innovations")
        residuals.append(residual / scale)

    h, l = residuals
    length = max(4, round(duration * fs))
    step = max(1, length // 4)
    candidates = int(round(t0 * fs)) + np.arange(-2, 3) * step
    controls = np.arange(fs, n - fs, length)
    controls = controls[np.abs(controls / fs - t0) > 1.5 + duration / 2]
    if len(controls) < 3:
        raise ValueError("insufficient off-source comparison windows")
    centers = np.r_[candidates, controls]
    starts = centers - length // 2
    ends = starts + length
    maxlag = int(np.floor(0.010 * fs))
    if starts.min() <= maxlag or ends.max() + maxlag >= n:
        raise ValueError("insufficient samples for requested comparison")
    eh = np.r_[0., np.cumsum(h*h)]
    el = np.r_[0., np.cumsum(l*l)]
    energy_h = eh[ends] - eh[starts]
    # Maximize |mean cross product| at physical delays. We also report the
    # normalized coherence and sign of that SAME extremum.
    best = np.full(len(centers), -np.inf)
    lags = np.zeros(len(centers), dtype=int)
    signs = np.ones(len(centers), dtype=int)
    coherences = np.zeros(len(centers))
    for lag in range(-maxlag, maxlag + 1):
        # Disjoint window indices ensure no circular wraparound for either series.
        lo = max(0, -lag)
        hi = min(n, n - lag)
        acc = np.r_[0., np.cumsum(h[lo:hi] * l[lo+lag:hi+lag])]
        dot = acc[ends-lo] - acc[starts-lo]
        val = np.abs(dot) / length
        improved = val > best
        if not improved.any():
            continue
        energy_l = el[ends+lag] - el[starts+lag]
        best[improved] = val[improved]
        coherences[improved] = np.clip(np.abs(dot[improved]) /
            np.sqrt(np.maximum(energy_h[improved] * energy_l[improved], 1e-30)), 0, 1)
        signs[improved] = np.sign(dot[improved])
        lags[improved] = lag
    winner = int(np.argmax(best[:len(candidates)]))
    background = np.median(best[len(candidates):])
    return Confirmation(first_second + left + centers[winner] / fs, lags[winner] / fs,
                        int(signs[winner]), float(best[winner]), float(coherences[winner]),
                        float(best[winner] / max(background, 1e-30)))


def sine_gaussian(fs, n_seconds, peak_time, frequency, q, hrss, phase=0.):
    """Detector-frame software injection with numerical integral h^2 dt = hrss^2.

    Envelope exp(-(pi*f*(t-t0)/Q)^2*2) uses sigma_t = Q/(2*pi*f).
    Both peak_time and the length are relative to the supplied crop.
    """
    if fs <= 0 or int(fs*n_seconds) != fs*n_seconds or not 0 < frequency < fs/2:
        raise ValueError("invalid injection sampling")
    if not 0 < peak_time < n_seconds or q <= 0 or hrss < 0:
        raise ValueError("invalid injection time, Q, or amplitude")
    t = np.arange(round(fs*n_seconds)) / fs - peak_time
    phase_time = 2 * np.pi * frequency * t
    wave = np.exp(-0.5 * (phase_time / q)**2) * np.cos(phase_time + phase)
    return wave * (hrss / np.sqrt(np.dot(wave, wave) / fs))
