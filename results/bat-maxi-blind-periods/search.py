"""Blind analysis of daily BAT and MAXI rate measurements. See METHODS.md."""
from functools import lru_cache
import hashlib

import numpy as np
from scipy import fft, linalg, optimize, signal, sparse

P_MIN = 2.05
P_MAX = 2000.
MC = 255
TRAIN = 64
ALPHA = .05
SURVEY_PASSES = 2
SYSTEMS = {
    "solar-day": (1., .00002), "sidereal-day": (.9972696, .00002),
    "Swift-precession": (53., 2.), "ISS-precession": (70., 2.),
    "synodic-Moon": (29.53059, .02), "sidereal-Moon": (27.32166, .02),
    "anomalistic-Moon": (27.55455, .02), "draconic-Moon": (27.21222, .02),
    "year": (365.256, .25),
}
INJ_PERIODS = (7.37, 41.81, 173.3, 683.3)
INJ_AMPLITUDES = (.10, .30, 1.)
INJ_PHASES = 8


def rng_for(source, purpose):
    seed = hashlib.sha256(f"bat-maxi-v1/{purpose}/{source}".encode()).digest()[:16]
    return np.random.default_rng(int.from_bytes(seed, "little"))


def survey_trials(n_grid, n_sources):
    return SURVEY_PASSES * int(n_grid) * int(n_sources)


def kish(errors):
    if not len(errors):
        return 0.
    w = (errors.min()/errors)**2
    return float(w.sum()**2/(w@w))


def bounded_uncertainties(t, e):
    """Limit design leverage using only observation times and declared errors."""
    if not len(e):
        return e, {"original_effective_samples":0., "effective_samples":0.,
                   "error_floor":None, "errors_raised":0, "weight_halves":[]}
    left = t <= (np.min(t)+np.max(t))/2
    groups = [e[left],e[~left]]
    floors = []
    for group in groups:
        if len(group):
            # At least ten equally precise observations (when present) and
            # at least one tenth of each half share the maximum possible weight.
            count = min(len(group),max(10,int(np.ceil(.10*len(group)))))
            floors.append(float(np.partition(group,count-1)[count-1]))
        else:
            floors.append(None)
    floor = max(x for x in floors if x is not None)
    corrected = np.maximum(e,floor)
    details = [{"samples":len(group),"order_statistic_error":f,
                "raw_kish":kish(group),"adjusted_kish":kish(np.maximum(group,floor))}
               for group,f in zip(groups,floors)]
    return corrected, {"original_effective_samples":kish(e),"effective_samples":kish(corrected),
                       "error_floor":floor,"errors_raised":int((e<floor).sum()),"weight_halves":details}


def prepare(curve):
    if curve is None:
        return None, {"status": "unavailable"}
    t, y, e = (np.asarray(curve, dtype=float) if isinstance(curve, tuple)
               else np.asarray(curve, dtype=float).T)
    keep = np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    t, y, e = t[keep], y[keep], e[keep]
    order = np.argsort(t)
    t, y, e = t[order], y[order], e[order]
    # Exact duplicate epochs cannot be counted as independent evidence.
    keep = np.r_[True, np.diff(t) > 0] if len(t) else np.zeros(0, bool)
    t, y, e = t[keep], y[keep], e[keep]
    e, weights = bounded_uncertainties(t,e)
    span = np.ptp(t) if len(t) else 0.
    meta = {"samples":len(t),"span_days":float(span),**weights}
    if len(t):
        meta.update(mjd_first=float(t[0]),mjd_last=float(t[-1]))
    # A curve can be scanned even if later two-half confirmation is impossible.
    if len(t)<16 or weights["effective_samples"]<6 or span<5*P_MIN:
        return None, dict(meta,status="insufficient_sampling")
    return (t, y, e), meta


class Continuum:
    """Nonnegative mixture of stationary OU processes plus independent errors."""

    def __init__(self, t, y, e, w):
        self.e, self.t = e, t
        span = np.ptp(t)
        self.timescales = np.unique(np.r_[0., np.geomspace(1., span/3, 8)])
        lags, empirical, uncertainty = [], [], []
        for lag in np.geomspace(1, span/3, 26):
            j = np.searchsorted(t, t + lag)
            j = np.clip(j, 1, len(t)-1)
            j -= (abs(t[j-1] - t - lag) < abs(t[j] - t - lag)).astype(int)
            dt = t[j] - t
            ok = (dt > 0) & (abs(dt-lag) <= max(.55, .15*lag))
            if ok.sum() < 15:
                continue
            v = (e[j[ok]]**2 + e[ok]**2) / 2
            pw = 1/v
            pw /= pw.sum()
            half_square = (y[j[ok]] - y[ok])**2 / 2
            lags.append(float(pw@dt[ok]))
            empirical.append(float(pw@(half_square-v)))
            # Prevent the sparsest, very long lags from dominating the fit.
            uncertainty.append(max(float(pw@half_square), float(pw@v))
                               / np.sqrt(min(100., 1/(pw@pw), span/lag)))
        scale = max(float(w@((y - w@y)**2)), float(w@(e**2)), np.finfo(float).tiny)
        if len(lags) >= 4:
            design = np.column_stack([np.ones(len(lags))] +
                [-np.expm1(-np.asarray(lags)/tau) for tau in self.timescales[1:]])
            weights = scale / np.maximum(uncertainty, scale*1e-10)
            solve = optimize.lsq_linear(design*weights[:,None], np.array(empirical)/scale*weights,
                                        bounds=(0, np.inf), method="bvls")
            self.variances = solve.x*scale
        else:
            self.variances = np.zeros_like(self.timescales)
            self.variances[0] = max(0., float(w@((y-w@y)**2-e**2)))
        self.enough = len(lags) >= 4
        self.method = "variogram"
        self.configure()

    def filters(self, tau):
        rho = np.exp(-np.diff(self.t)/tau)
        forward = np.zeros((2, len(self.t)))
        forward[0] = 1.
        forward[1, :-1] = -rho
        backward = np.zeros_like(forward)
        backward[1] = 1.
        backward[0, 1:] = -rho
        return forward, backward

    def configure(self):
        self.components = []
        self.white = self.e**2 + self.variances[0]
        scale = max(float(np.mean(self.white)), float(self.variances.sum()))
        for tau, var in zip(self.timescales[1:], self.variances[1:]):
            if var <= scale*1e-12:
                continue
            forward, backward = self.filters(tau)
            innovation = np.sqrt(np.r_[1., -np.expm1(-2*np.diff(self.t)/tau)] * var)
            self.components.append((var, forward, backward, innovation))

    def fit_spectrum(self, w, frequencies, powers, regions):
        """Fit median spectral levels using the native-time, window-folded model response."""
        regions = [r for r in regions if r.sum() >= 8]
        if len(regions) < 4:
            return
        # At seven representative frequencies per band orthonormalize the two
        # floating-mean sinusoid columns. The expected GLS power is then the
        # trace of their projected covariance. A forward recurrence suffices
        # for each quadratic form of an OU kernel on the irregular time grid.
        ids = np.array([np.linspace(np.flatnonzero(r)[0], np.flatnonzero(r)[-1],
                                   7).round().astype(int) for r in regions])
        angles = 2*np.pi*self.t[:,None]*frequencies[ids.ravel()]
        x = np.stack((np.cos(angles), np.sin(angles)), axis=-1)
        x -= np.einsum('n,nfi->fi', w, x, optimize=True)
        gram = np.einsum('n,nfi,nfj->fij', w, x, x, optimize=True)
        vals, vec = np.linalg.eigh(gram)
        basis = np.einsum('nfi,fij->nfj', x, vec, optimize=True)
        basis *= w[:,None,None] / np.sqrt(np.maximum(vals, 1e-12))[None,:,:]
        v = basis.reshape(len(w), -1)
        def band(values):
            return np.median(values.reshape(len(regions), 7, 2).sum(axis=2), axis=1)
        diagonal = np.sum(v*v, axis=0)
        measurement = band(np.sum((self.e[:,None]*v)**2, axis=0))
        response = [band(diagonal)]
        for tau in self.timescales[1:]:
            forward, _ = self.filters(tau)
            left = linalg.solve_banded((1,0), forward, v, check_finite=False)
            response.append(band(np.maximum(0., 2*np.sum(v*left,axis=0) - diagonal)))
        response = np.asarray(response).T
        observed = np.array([np.median(powers[r]) for r in regions])/np.log(2.)
        # The exponential median-to-mean conversion is approximate for red,
        # anisotropic projections. Independent forward simulations below check
        # the actual median bandpower distributions after fitting.
        width = np.array([r.sum() for r in regions])
        weights = np.sqrt(np.minimum(30., width/4)) / np.maximum.reduce([
            observed,measurement,np.full(len(observed),np.finfo(float).tiny)])
        scale = max(float(np.median(self.e)**2), float(np.median(observed/response[:,0])),
                    np.finfo(float).tiny)
        try:
            variance, _ = optimize.nnls(weights[:,None]*response*scale,
                weights*(observed-measurement), maxiter=1000)
        except (RuntimeError, ValueError, linalg.LinAlgError):
            return
        self.variances = variance*scale
        self.method = "window-folded spectral-band medians"
        self.enough = True
        self.configure()

    def draws(self, rng, count):
        vals = rng.standard_normal((len(self.e), count)) * np.sqrt(self.white)[:, None]
        for _, forward, _, innovation in self.components:
            x = rng.standard_normal(vals.shape) * innovation[:, None]
            vals += linalg.solve_banded((1, 0), forward, x, check_finite=False)
        return vals

    def bilinear(self, v):
        """v.T @ Cov(y) @ v for the actual times, including measurement errors."""
        cv = self.white[:, None]*v
        for var, forward, backward, _ in self.components:
            left = linalg.solve_banded((1, 0), forward, v, check_finite=False)
            right = linalg.solve_banded((0, 1), backward, v, check_finite=False)
            cv += var*(left+right-v)
        return v.T@cv


class Gridding:
    """Approximate nonuniform Fourier sums with a Gaussian interpolation kernel."""

    def __init__(self, elapsed, step, nfft):
        x = elapsed/step
        center = np.rint(x).astype(int)
        nodes = center[:, None] + np.arange(-5, 6)
        coefficients = np.exp(-.5*(nodes-x[:, None])**2)/np.sqrt(2*np.pi)
        self.spread = sparse.csr_matrix(
            (coefficients.ravel(), (nodes.ravel(), np.repeat(np.arange(len(x)), 11))),
            shape=(int(nodes.max())+1, len(x)))
        self.nfft, self.step = nfft, step

    def sums(self, values, indices):
        gridded = self.spread@values
        transform = fft.rfft(gridded.T, n=self.nfft, axis=-1, workers=1)
        deconvolve = np.exp(2*(np.pi*indices/self.nfft)**2)
        return transform[..., indices]*deconvolve


def extrema(scores):
    return signal.find_peaks(np.r_[-np.inf, scores, -np.inf])[0] - 1


def mc_p(maxima, value):
    return float((1 + np.count_nonzero(maxima >= value))/(1 + len(maxima)))


def design(t, frequency):
    angle = 2*np.pi*frequency*t
    return np.stack((np.ones_like(t), np.cos(angle), np.sin(angle)), axis=-1)


def fit(t, y, w, frequency):
    x = design(t, frequency)
    normal = x.T@(w[:, None]*x)
    b = linalg.pinvh(normal, check_finite=False) @ (x.T@(w*y))
    improvement = float(w@((x[:, 1:]@b[1:] - w@(x[:, 1:]@b[1:])/w.sum())**2))
    return b, improvement


class Search:
    def __init__(self, source_id, band, data, metadata, n_sources):
        t, self.y, e = data
        self.origin = float(t[0]-4)
        self.t, self.e = t-self.origin, e
        self.w = (e.min()/e)**2
        self.w /= self.w.sum()
        self.span = float(np.ptp(t))
        self.n_sources = n_sources
        self.scale = float(np.sqrt(self.w@((self.y - self.w@self.y)**2)))
        self.noise = Continuum(self.t, self.y, e, self.w)
        nfft = fft.next_fast_len(int(np.ceil(8*(self.span+10))), real=True)
        self.df = 2/nfft
        self.idx = np.arange(int(np.ceil(max(4/self.span, 1/P_MAX)/self.df)),
                             int(np.floor((1/P_MIN)/self.df))+1)
        self.freq = self.idx*self.df
        self.grid = Gridding(self.t, .5, nfft)
        fine = Gridding(self.t, .25, 2*nfft)
        s1 = fine.sums(self.w, self.idx)
        s2 = fine.sums(self.w, 2*self.idx)
        cc = (1+s2.real)/2 - s1.real**2
        ss = (1-s2.real)/2 - s1.imag**2
        cs = s2.imag/2 - s1.real*s1.imag
        det = cc*ss-cs**2
        good = det > 1e-8
        self.gcc = np.divide(ss, det, where=good, out=np.zeros_like(ss))
        self.gss = np.divide(cc, det, where=good, out=np.zeros_like(cc))
        self.gcs = np.divide(-2*cs, det, where=good, out=np.zeros_like(cs))

        win_idx = np.arange(1, int(1.45/self.df))
        ex = np.abs(fine.sums(self.w, win_idx))**2
        sampling = np.abs(fine.sums(np.full(len(t), 1/len(t)), win_idx))**2
        win = np.maximum(ex, sampling)
        win_peaks, _ = signal.find_peaks(win, height=.01, distance=max(1, int(.7/self.span/self.df)))
        # Include the finite-baseline sidelobes just beyond one resolution element.
        win_peaks = win_peaks[win_idx[win_peaks]*self.df > 1.1/self.span]
        self.measured_window_f = win_idx[win_peaks]*self.df
        self.measured_window_p = win[win_peaks]
        self.spacecraft = "Swift-precession" if band=="bat" else "ISS-precession"
        other_spacecraft = "ISS-precession" if band=="bat" else "Swift-precession"
        all_windows = [(k,1/p,dp/p**2) for k,(p,dp) in SYSTEMS.items()]
        self.windows = [item for item in all_windows if item[0]!=other_spacecraft]
        self.advisories = []
        for label,fw,width in all_windows:
            if label==other_spacecraft:
                for h in range(1,5):
                    self.advisories.append((f"{h} x {label} (other spacecraft)",h*fw,h*width))
        for i,(label,f1,d1) in enumerate(all_windows):
            for other,f2,d2 in all_windows[i+1:]:
                for side in (abs(f1-f2),f1+f2):
                    if self.freq[0]-1/self.span <= side <= self.freq[-1]+1/self.span:
                        self.advisories.append((label+" x "+other,side,d1+d2))
        self.windows += [("observed-window",float(fw),0.) for fw in self.measured_window_f]
        self.direct = np.zeros(len(self.freq),bool)
        for label,fw,width in self.windows:
            # A theoretical combination by itself is advisory. Any actual
            # window line in its tolerance range is vetoed at its measured
            # frequency, not throughout the nominal range of the combination.
            for h in (range(1,5) if label!="observed-window" else (1,)):
                self.direct |= np.abs(self.freq-h*fw) < 1/self.span + h*width

        edges = np.geomspace(self.freq[0], self.freq[-1]*(1+1e-8), 12)
        regions = [(self.freq >= low)&(self.freq < high) for low, high in zip(edges[:-1], edges[1:])]
        regions = [r for r in regions if r.sum() >= 8]
        if not regions:
            regions = [np.ones(len(self.freq),dtype=bool)]
        self.base = self.transforms(self.y)[0]
        self.noise.fit_spectrum(self.w, self.freq, self.power(self.base), regions)
        # Use independent draws for fitting the score normalization and its null distribution.
        train_rng = rng_for(source_id, band+"-train")
        self.level = np.zeros(len(self.freq))
        for start in range(0, TRAIN, 8):
            p = self.raw(self.noise.draws(train_rng, min(8, TRAIN-start)))
            self.level += p.sum(axis=0)
        self.level = np.maximum(self.level/TRAIN, np.finfo(float).tiny)
        self.observed = self.scores(self.base)
        self.maxima, bands = [], []
        null_rng = rng_for(source_id, band+"-calibrate")
        for start in range(0, MC, 8):
            p = self.raw(self.noise.draws(null_rng, min(8, MC-start))) / self.level
            self.maxima.extend(p.max(axis=1))
            bands.extend(np.column_stack([np.median(p[:, r], axis=1) for r in regions]))
        self.maxima = np.asarray(self.maxima)
        low, high = np.quantile(bands, [.01, .99], axis=0)
        observed_bands = np.array([np.median(self.observed[r]) for r in regions])
        bad = (observed_bands < .5*low) | (observed_bands > 2*high)
        self.model_ok = bool(self.noise.enough and len(regions)>=2 and bad.sum()<2)
        self.alias_offsets = []
        for label, fw, width in self.windows:
            for h in (range(1, 5) if label != "observed-window" else (1,)):
                self.alias_offsets.append((label, h*fw, 1/self.span+h*width))
        self.alias_frequencies = np.array([item[1] for item in self.alias_offsets])
        self.alias_widths = np.array([item[2] for item in self.alias_offsets])
        self.threshold = np.sort(self.maxima)[-2]  # rank p < .01 with 255 null draws
        self.half_masks = [self.t <= (self.t[0]+self.t[-1])/2, self.t > (self.t[0]+self.t[-1])/2]
        self.meta = dict(metadata, status="searched", min_period_days=float(1/self.freq[-1]),
                         max_period_days=float(1/self.freq[0]), frequencies=len(self.freq),
                         grid_spacing_per_day=float(self.df), weighted_rms=self.scale,
                         sampling_window_lines=len(win_peaks), model_ok=self.model_ok,
                         survey_passes=SURVEY_PASSES, survey_trial_factor=survey_trials(len(self.freq),n_sources),
                         native_precession=self.spacecraft,
                         red_fit_method=self.noise.method, red_timescale_days=self.noise.timescales.tolist(),
                         red_variances=self.noise.variances.tolist(),
                         continuum_diagnostics=[{"period_days": [float(1/self.freq[r][-1]),float(1/self.freq[r][0])],
                             "observed_median": float(a), "null_01": float(b), "null_99": float(c),
                             "outside_tolerance": bool(d)} for r,a,b,c,d in zip(regions,observed_bands,low,high,bad)],
                         training_simulations=TRAIN, calibration_simulations=MC)

    def transforms(self, y):
        y = np.asarray(y).reshape(len(self.t), -1)
        return self.grid.sums(self.w[:, None]*(y-(self.w@y)[None, :]), self.idx)

    def power(self, sums):
        a, b = sums.real, sums.imag
        return np.maximum(0., self.gcc*a*a+self.gss*b*b+self.gcs*a*b)

    def raw(self, y):
        return self.power(self.transforms(y))

    def scores(self, sums):
        return self.power(sums)/self.level

    @lru_cache(maxsize=96)
    def precision(self, k):
        x = design(self.t, self.freq[k])[:, 1:]
        v = self.w[:, None]*(x-self.w@x)
        cov = self.noise.bilinear(v)
        return v, linalg.pinvh(cov, check_finite=False)

    def analytic(self, k, y):
        v, inv = self.precision(int(k))
        z = v.T@y
        logp = float(-.5*(z@inv@z))
        return logp, min(1., float(np.exp(max(-745., logp+np.log(survey_trials(len(self.idx),self.n_sources))))))

    def clock_advisories(self,k):
        frequency = self.freq[k]
        near = np.abs(self.measured_window_f-frequency) < 1/self.span
        peak = float(np.max(self.measured_window_p[near])) if near.any() else None
        return [{"clock":label,"predicted_frequency":float(fw),"measured_window_power":peak}
                for label,fw,width in self.advisories if abs(frequency-fw)<1/self.span+width]

    def window_reason(self, k, parents):
        f, res = self.freq[k], 1/self.span
        if self.direct[k]:
            return "direct instrumental or sampling window"
        for other in parents:
            if k == other:
                continue
            g = self.freq[other]
            if abs(f-g) < 1.25*res:
                return "unresolved stronger peak"
            for h in range(2, 7):
                if abs(f-h*g) < 1.25*res:
                    return f"harmonic {h} of stronger peak"
            offset = self.alias_frequencies
            found = (np.minimum(np.abs(f-np.abs(g-offset)), np.abs(f-(g+offset))) < self.alias_widths)
            if np.any(found):
                label = self.alias_offsets[int(np.argmax(found))][0]
                return f"sideband of stronger peak at {1/g:.6g} d ({label})"
        return None

    def ranked(self, scores):
        local = extrema(scores)
        selected = local[scores[local] > self.threshold]
        selected = selected[np.argsort(-scores[selected], kind="stable")]
        blocked = np.zeros(len(self.freq), dtype=int)
        rows = []
        for k in selected:
            if self.direct[k]:
                reason = "direct native clock/harmonic or measured sampling line"
            elif blocked[k]:
                reason = self.window_reason(k, [blocked[k]-1]) or "sideband of stronger peak"
            else:
                reason = None
            rows.append((int(k), reason))
            if reason and not self.direct[k]:
                continue
            # Mark intervals once per parent rather than comparing all pairs of maxima.
            g = self.freq[k]
            centers = np.r_[g, np.arange(2,7)*g, np.abs(g-self.alias_frequencies), g+self.alias_frequencies]
            widths = np.r_[np.full(6,1.25/self.span), self.alias_widths, self.alias_widths]
            lower = np.searchsorted(self.freq, centers-widths, side="right")
            upper = np.searchsorted(self.freq, centers+widths, side="left")
            delta = (np.bincount(lower, minlength=len(self.freq)+1) -
                     np.bincount(upper, minlength=len(self.freq)+1))
            affected = (np.cumsum(delta[:-1]) > 0) & (blocked == 0)
            affected[k] = False
            blocked[affected] = k+1
        below = local[scores[local] <= self.threshold]
        vetoed = int(np.count_nonzero(self.direct[below] | (blocked[below] > 0)))
        counts = {"local_maxima": len(local), "direct_window_maxima": int(self.direct[local].sum()),
                  "window_or_alias_maxima": vetoed + sum(reason is not None for _,reason in rows),
                  "above_mc_threshold": len(selected)}
        return rows, counts

    def peaks(self, sums, y):
        scores = self.scores(sums)
        ranks, counts = self.ranked(scores)
        rows = []
        for k, reason in ranks:
            logp, bound = self.analytic(k, y)
            rows.append({"index": k, "period_days": float(1/self.freq[k]), "score": float(scores[k]),
                "source_global_mc_fap": mc_p(self.maxima, scores[k]),
                "single_grid_log_p": logp, "survey_fap_bound": bound,
                "window_reason": reason, "model_ok": self.model_ok,
                "advisory_clock_proximity": self.clock_advisories(k)})
        return rows, counts

    def halves(self, frequency, y):
        result = []
        for mask in self.half_masks:
            t, yy, w = self.t[mask], y[mask], self.w[mask]
            if len(t) < 16 or np.ptp(t)*frequency < 1.9 or w.sum()**2/(w@w) < 10:
                return {"pass": False, "reason": "insufficient half coverage", "halves": result}
            w = w/w.sum()
            span = np.ptp(t)
            # Search a wider interval than the acceptance interval.
            local = frequency + np.arange(-12,13)/(6*span)
            local = local[(local > 0) & (local < .5)]
            b, _ = fit(t, yy, w, frequency)
            peak = max(local, key=lambda f: fit(t, yy, w, f)[1])
            phase = float((np.arctan2(b[2], b[1])/(2*np.pi)) % 1)
            # In-band/half detection uses the same red covariance, at the full-data period.
            x = design(t, frequency)[:, 1:]
            vv = np.zeros((len(self.t), 2))
            vv[mask] = w[:, None]*(x-w@x)
            z = vv[mask].T@yy
            p = float(np.exp(-.5*(z@linalg.pinvh(self.noise.bilinear(vv))@z)))
            result.append({"period_days": float(1/peak), "frequency_offset_rayleigh": float(abs(peak-frequency)*span),
                           "semiamplitude": float(np.hypot(*b[1:])), "phase_cycles": phase,
                           "local_red_p": p})
        drift = abs((result[0]["phase_cycles"] - result[1]["phase_cycles"] + .5) % 1 - .5)
        ok = all(h["frequency_offset_rayleigh"] <= 1.01 and h["local_red_p"] <= .05 for h in result)
        return {"pass": bool(ok and drift <= .20), "phase_drift_cycles": drift, "halves": result}

    def describe(self, row):
        k = row["index"]
        f0 = self.freq[k]
        lo, hi = max(self.freq[0], f0-self.df), min(self.freq[-1], f0+self.df)
        f = optimize.minimize_scalar(lambda f: -fit(self.t, self.y, self.w, f)[1],
                                     bounds=(lo,hi), method="bounded",
                                     options={"xatol": self.df/10000}).x
        b, _ = fit(self.t, self.y, self.w, f)
        amp = np.hypot(*b[1:])
        phase = 2*np.pi*self.t*f
        x = np.column_stack((design(self.t, f), 2*np.pi*self.t*(-b[1]*np.sin(phase)+b[2]*np.cos(phase))))
        v = self.w[:, None]*x
        # Scale design columns before inversion for rate units ranging over many orders of magnitude.
        scale = np.sqrt(np.maximum(np.sum(v*x, axis=0), np.finfo(float).tiny))
        v = v/scale
        inv = linalg.pinvh(v.T@(x/scale))
        cov = (inv@self.noise.bilinear(v)@inv) / np.outer(scale, scale)
        grad = b[1:]/max(amp, np.finfo(float).tiny)
        folds = []
        bin_ids = np.minimum(15, ((self.t*f % 1)*16).astype(int))
        for i in range(16):
            m = bin_ids == i
            if not m.any():
                folds.append({"phase": (i+.5)/16, "n": 0})
                continue
            w = self.w[m]/self.w[m].sum()
            folds.append({"phase": (i+.5)/16, "n": int(m.sum()), "rate": float(w@self.y[m]),
                          "measurement_se": float(np.sqrt((w**2)@(self.e[m]**2)))})
        return dict(row, period_days=float(1/f), uncertainty_days=float(np.sqrt(max(0.,cov[3,3]))/f**2),
                    semiamplitude=float(amp), semiamplitude_se=float(np.sqrt(max(0.,grad@cov[1:3,1:3]@grad))),
                    phase_origin_mjd=self.origin, fold=folds)


def accepted(row):
    return row["window_reason"] is None and row["survey_fap_bound"] <= ALPHA and row["model_ok"]


def match(bat, maxi, br, mr):
    # Half a Rayleigh bin from each independent full-baseline search.
    return abs(1/br["period_days"]-1/mr["period_days"]) <= .5/bat.span + .5/maxi.span


def inject(source, instruments):
    # Stratified by stable hash, never by the strength or absence of an observed peak.
    rng = rng_for(source, "injections")
    assignments = [p for p in INJ_PERIODS if all(s.freq[0] <= 1/p <= s.freq[-1] for s in instruments.values())]
    if not assignments:
        return []
    period = float(rng.choice(assignments))
    phases = rng.uniform(0, 1, INJ_PHASES)
    shared_origin = min(s.origin for s in instruments.values())
    prepared = {}
    for name, s in instruments.items():
        angle = 2*np.pi*((s.origin-shared_origin+s.t)/period)
        sincos = np.column_stack((np.sin(angle),np.cos(angle)))
        prepared[name] = (sincos, s.transforms(sincos))
    rows = []
    for amplitude in INJ_AMPLITUDES:
        hits = {key: 0 for key in instruments}
        both = 0
        for phase in phases:
            found = {}
            for name, s in instruments.items():
                dx, ds = prepared[name]
                phi = 2*np.pi*phase
                b = amplitude*s.scale*np.array([np.cos(phi),np.sin(phi)])
                y = s.y + dx@b
                scores = s.scores(s.base + b@ds)
                nearby = np.flatnonzero(np.abs(s.freq-1/period) <= 1/s.span)
                if not s.model_ok or not len(nearby) or scores[nearby].max() <= s.threshold:
                    continue
                for k, window in s.ranked(scores)[0]:
                    if window is not None or abs(s.freq[k]-1/period) > 1/s.span:
                        continue
                    _, bound = s.analytic(k, y)
                    if bound <= ALPHA and s.halves(s.freq[k], y)["pass"]:
                        found.setdefault(name, []).append({"period_days": float(1/s.freq[k])})
                hits[name] += name in found
            if "bat" in found and "maxi" in found:
                both += any(match(instruments["bat"], instruments["maxi"], b, m)
                            for b in found["bat"] for m in found["maxi"])
        rows.append({"period_days": period, "relative_semiamplitude": amplitude, "trials": len(phases),
                     "rate_semiamplitudes": {name: amplitude*s.scale for name,s in instruments.items()},
                     "recovered_by_band": hits, "joint_recovered": both if len(instruments)==2 else None})
    return rows


def analyze(source_id, bat, maxi=None, n_sources=2000):
    """Return JSON-compatible results for one source, using full-roster multiplicity."""
    out = {"source_id": str(source_id), "instruments": {}, "candidates": [], "injections": [], "pairs": []}
    instruments, peaks = {}, {}
    for name, curve in (("bat", bat), ("maxi", maxi)):
        cleaned, meta = prepare(curve)
        out["instruments"][name] = meta
        if cleaned is None:
            continue
        s = Search(source_id, name, cleaned, meta, n_sources)
        instruments[name] = s
        rows, counts = s.peaks(s.base, s.y)
        out["instruments"][name] = dict(s.meta, **counts, peaks=rows)
        peaks[name] = [r for r in rows if accepted(r)]
    if not instruments:
        out["status"] = "no_searchable_curves"
        return out
    if len(instruments) == 2:
        bs, ms = instruments["bat"], instruments["maxi"]
        for b in peaks["bat"]:
            for m in peaks["maxi"]:
                if not match(bs, ms, b, m):
                    continue
                bh, mh = bs.halves(bs.freq[b["index"]],bs.y), ms.halves(ms.freq[m["index"]],ms.y)
                pair = {"bat_grid_period_days": b["period_days"], "maxi_grid_period_days": m["period_days"],
                        "bat_stability": bh, "maxi_stability": mh, "pass": bool(bh["pass"] and mh["pass"])}
                out["pairs"].append(pair)
                if pair["pass"]:
                    bd, md = bs.describe(b), ms.describe(m)
                    out["candidates"].append({"period_days": bd["period_days"], "uncertainty_days": bd["uncertainty_days"],
                        "bat": bd, "maxi": md, "stability": pair,
                        "survey_fap_bound": max(bd["survey_fap_bound"], md["survey_fap_bound"])})
    out["injections"] = inject(source_id, instruments)
    out["status"] = "dual_search" if len(instruments)==2 else "single_instrument_search"
    # lru_cache would otherwise retain Search instances between sources in the shard.
    Search.precision.cache_clear()
    return out
