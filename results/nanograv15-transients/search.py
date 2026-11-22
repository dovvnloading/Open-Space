"""Predeclared narrowband single-pulsar transient search. See SINGLE_SEARCH.md."""
import argparse
import hashlib
import json
import traceback
from pathlib import Path

import numpy as np
from scipy import linalg, optimize, stats

YEAR = 365.25
SINGLE_ALPHA = 0.005
PULSARS = 68
MAX_EPOCHS = 256
INJECTIONS = 256
WIDTHS = (14., 60., 240.)
FIXED = ("phase_step", "spin_step_or_memory", "spin_slope_step", "dm_step")
SHAPES = ("achromatic_pulse", "shapiro_log", "doppler_0", "doppler_45",
          "doppler_90", "doppler_135", "dm_pulse", "dm_tail", "profile_like_tail")


def basis(matrix):
    norms = np.linalg.norm(matrix, axis=0)
    matrix = matrix[:, norms > 0] / norms[norms > 0]
    if not matrix.shape[1]:
        return matrix
    q, r, _ = linalg.qr(matrix, mode="economic", pivoting=True, check_finite=False)
    diag = np.abs(np.diag(r))
    rank = np.sum(diag > max(matrix.shape) * np.finfo(float).eps * diag.max())
    return q[:, :rank]


def remove(q, x):
    return x - q @ (q.T @ x)


class Projection:
    """Inner products after ECORR whitening and timing/noise nuisance projection."""

    def __init__(self, data, fit_soft):
        self.t = np.asarray(data["t"], dtype=float)
        self.y = np.asarray(data["y"], dtype=float)
        self.freq = np.asarray(data["freq"], dtype=float)
        self.err = np.asarray(data["err"], dtype=float)
        self.backend = np.asarray(data["backend"], dtype=str)
        self.ra = float(data["ra"])
        self.dec = float(data["dec"])
        n = len(self.t)
        if n < 16 or np.ptp(self.t) <= 0:
            raise ValueError("insufficient time coverage")
        if np.any(self.freq <= 0) or np.any(self.err <= 0):
            raise ValueError("nonpositive error or frequency")
        design = np.asarray(data["design"], dtype=float).reshape(n, -1)
        if not all(np.isfinite(x).all() for x in (self.t, self.y, self.freq, self.err, design)):
            raise ValueError("nonfinite reduction")
        self.chrom = (1400 / self.freq) ** 2
        ids = np.asarray(data.get("ecorr_group", np.full(n, -1)))
        ecorr = np.asarray(data.get("ecorr_us", np.zeros(n)), dtype=float)
        if np.any(ecorr < 0) or not np.isfinite(ecorr).all():
            raise ValueError("invalid ECORR noise")
        if np.any((ids < 0) & (ecorr > 0)):
            raise ValueError("nonzero ECORR without group")
        self.blocks = []
        for group in np.unique(ids[ids >= 0]):
            ix = np.flatnonzero(ids == group)
            v = ecorr[ix] / self.err[ix]
            v2 = v @ v
            if v2:
                self.blocks.append((ix, v, -1 / (np.sqrt(1 + v2) * (1 + np.sqrt(1 + v2)))))
        backend_design = self.backend[:, None] == np.unique(self.backend)
        x = (self.t - self.t.mean()) / np.ptp(self.t)
        self.q_timing = basis(self.weight(np.column_stack(
            (design, backend_design, x, x*x, self.chrom, self.chrom*x))))
        budget = n - self.q_timing.shape[1] - 14
        self.modes = min(20, max(1, int(np.ptp(self.t) // 180)), max(0, budget // 4))
        h = np.arange(1, self.modes + 1)
        phase = 2 * np.pi * x[:, None] * h
        modes = np.column_stack((np.sin(phase), np.cos(phase)))
        f = self.weight(np.column_stack((modes, self.chrom[:, None] * modes)))
        self.f = remove(self.q_timing, f)
        q_fourier = basis(self.f)
        # A second projection prevents tiny nuisance overlaps in rank-revealing QR.
        if q_fourier.shape[1]:
            singular = linalg.svd(remove(self.q_timing, q_fourier), full_matrices=False)
            q_fourier = singular[0][:, singular[1] > 1.e-7]
        self.q = np.column_stack((self.q_timing, q_fourier))
        self.rank = self.q.shape[1]
        self.dof = n - self.rank
        if self.dof < 9:
            raise ValueError("fewer than nine degrees of freedom after nuisance projection")
        self.ry = self.whiten(self.y)
        self.rss = float(self.ry @ self.ry)
        if self.rss <= 0:
            raise ValueError("no measurable postfit residual variance")
        self.white_scale = np.sqrt(self.rss / self.dof)
        self.noise = {"scale": float(self.white_scale), "fourier_modes_per_process": self.modes,
                      "guard_rank": self.rank, "guard_dof": self.dof,
                      "ecorr_blocks": len(self.blocks)}
        self.dm_restored = bool(data.get("dm_restored", False))
        self.dm_retained = bool(data.get("dm_retained", self.dm_restored))
        self.dm = np.asarray(data.get("dm", np.empty(0)), dtype=float)
        self.dm_err = np.asarray(data.get("dm_err", np.empty(0)), dtype=float)
        self.u_soft = np.empty((n, 0))
        self.soft_delta = np.zeros(0)
        if fit_soft:
            self.fit_processes()

    def weight(self, x):
        w = x / (self.err if x.ndim == 1 else self.err[:, None])
        # (I + v v^T)^(-1/2) from its sole nontrivial eigenvalue.
        for ix, v, factor in self.blocks:
            inner = v @ w[ix]
            if x.ndim == 1:
                w[ix] += factor * v * inner
            else:
                w[ix] += factor * v[:, None] * inner
        return w

    def whiten(self, x):
        return remove(self.q, self.weight(np.asarray(x)))

    def soft_whiten(self, x):
        w = remove(self.q_timing, self.weight(np.asarray(x)))
        coeff = self.u_soft.T @ w
        multiplier = self.soft_delta if w.ndim == 1 else self.soft_delta[:, None]
        return w + self.u_soft @ (multiplier * coeff)

    def fit_processes(self):
        if not self.modes:
            return
        f = self.f
        z = remove(self.q_timing, self.weight(self.y))
        zz = z @ z
        gram = f.T @ f
        rhs = f.T @ z
        dof = len(z) - self.q_timing.shape[1]
        log_k = np.log(np.tile(np.arange(1, self.modes + 1), 2))
        identity = np.eye(len(rhs))
        def amplitudes(theta):
            return np.exp(np.r_[theta[0] - .5 * theta[2] * log_k,
                                theta[1] - .5 * theta[3] * log_k])
        def objective(theta, details=False):
            a = amplitudes(theta)
            chol = linalg.cho_factor(identity + a[:, None] * gram * a, check_finite=False)
            vz = a * rhs
            ss = max(float(zz - vz @ linalg.cho_solve(chol, vz, check_finite=False)),
                     np.finfo(float).tiny)
            if details:
                return np.sqrt(ss / dof)
            return dof * np.log(ss / dof) + 2 * np.log(np.diag(chol[0])).sum()
        # Bound amplitudes in microseconds per unit white-noise scale.
        start = np.clip(np.log(np.median(self.err)), -8, 8)
        opt = optimize.minimize(objective, [start, start, 4., 8 / 3], method="L-BFGS-B",
                                bounds=[(-12., 12.), (-12., 12.), (0., 7.), (0., 7.)],
                                options={"maxiter": 75, "ftol": 1.e-9})
        a = amplitudes(opt.x)
        self.u_soft, s, _ = linalg.svd(f * a, full_matrices=False, check_finite=False)
        self.soft_delta = 1 / np.hypot(1, s) - 1
        scale = objective(opt.x, details=True)
        self.noise.update({"reml_converged": bool(opt.success), "reml_scale": float(scale),
                           "red_fundamental_coefficient_us": float(scale * np.exp(opt.x[0])),
                           "dm_fundamental_coefficient_us_at_1400": float(scale * np.exp(opt.x[1])),
                           "red_power_index": float(opt.x[2]), "dm_power_index": float(opt.x[3])})


def fit_pulsar(path, fit_soft=True):
    with np.load(path, allow_pickle=False) as arrays:
        return Projection(arrays, fit_soft)


def waveform(p, name, epoch, width):
    """Return waveform for unit amplitude in microseconds. Columns are epochs."""
    d = p.t[:, None] - np.atleast_1d(epoch)
    after = d >= 0
    chrom = p.chrom[:, None]
    if name == "phase_step":
        return after.astype(float)
    if name == "dm_step":
        return chrom * after
    if name == "spin_step_or_memory":
        return np.maximum(d, 0) / YEAR
    if name == "spin_slope_step":
        return (np.maximum(d, 0) / YEAR) ** 2
    x = d / width
    if name == "shapiro_log":
        return -np.log1p(x*x)
    if name.startswith("doppler_"):
        angle = np.deg2rad(float(name.split("_")[1]))
        return np.cos(angle) * (np.hypot(1, x) - 1) - np.sin(angle) * np.arcsinh(x)
    if name == "achromatic_pulse":
        return np.exp(-.5*x*x)
    if name == "dm_pulse":
        return np.exp(-.5*x*x) * chrom
    tail = np.exp(-np.maximum(x, 0)) * after
    if name == "dm_tail":
        return tail * chrom
    if name == "profile_like_tail":
        return tail * chrom**2
    raise ValueError(name)


def epoch_grid(t):
    dates = np.unique(np.floor(t))
    if len(dates) < 12:
        raise ValueError("fewer than 12 distinct observing days")
    mids = (dates[:-1] + dates[1:] + 1) / 2
    mids = mids[3:-3]
    indices = np.unique(np.linspace(0, len(mids)-1, min(MAX_EPOCHS, len(mids))).round().astype(int))
    return mids[indices]


def score(p, residualized):
    """Two-sided exact t-test given the fixed nuisance basis and relative covariance."""
    norm = np.linalg.norm(residualized, axis=0)
    unit = residualized / np.where(norm > 0, norm, 1)
    dot = p.ry @ unit
    variance = np.maximum(p.rss - dot*dot, np.finfo(float).tiny) / (p.dof - 1)
    t = dot / np.sqrt(variance)
    amp = dot / np.where(norm > 0, norm, 1)
    se = np.sqrt(variance) / np.where(norm > 0, norm, 1)
    return norm, unit, t, amp, se


def measurable(p, raw, norm):
    # Relative cutoff suppresses numerically unidentifiable timing degeneracies.
    return norm > np.maximum(1.e-14, 1.e-7 * np.linalg.norm(p.weight(raw), axis=0))


def fixed_check(p, event):
    epoch = event["epoch_mjd"]
    dates = np.unique(np.floor(p.t))
    if np.sum(dates < epoch) < 4 or np.sum(dates > epoch) < 4:
        return {"status": "epoch_outside_supported_span"}
    raw = waveform(p, event["template"], epoch, event["width_d"])
    norm, _, t, amp, se = score(p, p.whiten(raw))
    if not measurable(p, raw, norm)[0]:
        return {"status": "degenerate"}
    same = bool(amp[0] * event["amplitude_us"] > 0)
    # Direction fixed by the independent discovery data for EPTA.
    signed_t = np.sign(event["amplitude_us"]) * t[0]
    return {"status": "tested", "same_sign": same, "amplitude_us": float(amp[0]),
            "se_us": float(se[0]), "t": float(t[0]),
            "directional_p": float(stats.t.sf(signed_t, p.dof - 1)),
            "dm_restored": p.dm_restored, "dm_retained": p.dm_retained}


def wideband_dm_check(p, event):
    if not p.dm_retained or p.dm.shape != p.t.shape or p.dm_err.shape != p.t.shape:
        return {"status": "direct_dm_unavailable"}
    ok = np.isfinite(p.dm) & np.isfinite(p.dm_err) & (p.dm_err > 0)
    days = np.floor(p.t)
    # Use one measurement per day to avoid counting same-visit DMs as replicates.
    ix = []
    for day in np.unique(days[ok]):
        selection = np.flatnonzero(ok & (days == day))
        ix.append(selection[np.argmin(p.dm_err[selection])])
    ix = np.asarray(ix, dtype=int)
    t = p.t[ix]
    epoch = event["epoch_mjd"]
    if np.sum(t < epoch) < 4 or np.sum(t > epoch) < 4:
        return {"status": "unsupported_epoch"}
    x = (t - t.mean()) / np.ptp(t)
    instruments = p.backend[ix, None] == np.unique(p.backend[ix])
    base = np.column_stack((x, x*x, instruments))
    modes = min(p.modes, max(0, (len(t) - base.shape[1] - 10)//2))
    phase = 2 * np.pi * x[:, None] * np.arange(1, modes+1)
    sigma = p.dm_err[ix]
    q = basis(np.column_stack((base, np.sin(phase), np.cos(phase))) / sigma[:, None])
    dof = len(t) - q.shape[1] - 1
    if dof < 4:
        return {"status": "insufficient_dm_degrees_of_freedom"}
    y = remove(q, p.dm[ix] / sigma)
    ach = waveform(p, event["template"], epoch, event["width_d"]).ravel()[ix] / p.chrom[ix]
    h = remove(q, ach/sigma)
    norm = np.linalg.norm(h)
    if norm <= 1.e-7 * np.linalg.norm(ach/sigma):
        return {"status": "degenerate"}
    a = float(h @ y / norm**2)
    scale = np.linalg.norm(y - a*h) / np.sqrt(dof)
    se = max(float(scale / norm), np.finfo(float).tiny)
    return {"status": "tested", "distinct_days": len(ix), "dof": dof,
            "delta_dm_pc_cm3": a, "se_pc_cm3": se, "scale": float(scale),
            "same_sign": bool(a*event["amplitude_us"] > 0),
            "directional_p": float(stats.t.sf(np.sign(event["amplitude_us"])*a/se, dof))}


def group_crossings(events):
    # Display neighborhoods never alter validation, significance or raw rows.
    unassigned = set(range(len(events)))
    groups = []
    for i in sorted(unassigned, key=lambda j: events[j]["p_local"]):
        if i not in unassigned:
            continue
        t0 = events[i]["epoch_mjd"]
        members = sorted(j for j in unassigned if abs(events[j]["epoch_mjd"] - t0) <= 30.)
        unassigned.difference_update(members)
        groups.append({"representative_index": i, "primary_crossing_indices": members,
                       "survivor_indices": [j for j in members if not events[j]["unpassed"]]})
    return groups


def chromatic_fit(p, raw, name):
    """Fit coincident copies of this shape with three different spectra."""
    expected = 1 if name.startswith("dm_") else (2 if name == "profile_like_tail" else 0)
    reference = raw / p.chrom**expected
    templates = p.whiten(reference[:, None] * p.chrom[:, None]**np.arange(3))
    norms = np.linalg.norm(templates, axis=0)
    if np.any(norms == 0):
        return {"status": "degenerate"}
    u, s, vt = linalg.svd(templates/norms, full_matrices=False, check_finite=False)
    if s[-1] < 1.e-6*s[0]:
        return {"status": "unresolved_chromaticity"}
    coef = vt.T @ ((u.T @ p.ry)/s)
    unexplained = p.ry - u @ (u.T @ p.ry)
    scale = np.linalg.norm(unexplained)/np.sqrt(p.dof-3)
    uncertainty = np.sqrt(np.sum((vt.T/s)**2, axis=1)) * scale
    partial_t = coef[expected]/max(uncertainty[expected], np.finfo(float).eps)
    return {"status": "tested", "joint_amplitudes_us_nu0_nu2_nu4": (coef/norms).tolist(),
            "joint_se_us_nu0_nu2_nu4": (uncertainty/norms).tolist(),
            "expected_exponent_partial_t": float(partial_t),
            "expected_exponent_partial_p": float(2*stats.t.sf(abs(partial_t), p.dof-3))}


def delete_visit(p, raw):
    """Delete the day with largest absolute contribution to the template score."""
    z = p.whiten(raw).ravel()
    day = np.floor(p.t).astype(int)
    days, inv = np.unique(day, return_inverse=True)
    influence = np.abs(np.bincount(inv, z * p.ry))
    chosen = days[np.argmax(influence)]
    mask = day == chosen
    # Preserve entire ECORR blocks across a midnight boundary.
    for ix, _, _ in p.blocks:
        if np.any(mask[ix]):
            mask[ix] = True
    q = p.q[mask]
    w = np.column_stack((p.ry[mask], z[mask]))
    matrix = np.eye(mask.sum()) - q @ q.T
    eig = linalg.eigvalsh(matrix)
    if eig[0] < 1.e-8 or p.dof-mask.sum() < 3:
        return {"status": "deletion_loses_nuisance_rank", "mjd_day": int(chosen)}
    delta = w.T @ linalg.solve(matrix, w, assume_a="pos", check_finite=False)
    zz = z @ z - delta[1, 1]
    yz = p.ry @ z - delta[0, 1]
    if zz <= 1.e-10 * (z @ z):
        return {"status": "single_visit_support", "mjd_day": int(chosen)}
    df = p.dof - int(mask.sum()) - 1
    variance = max(p.rss - delta[0, 0] - yz**2/zz, np.finfo(float).tiny) / df
    t = yz / np.sqrt(variance*zz)
    return {"status": "tested", "mjd_day": int(chosen), "amplitude_us": float(yz/zz),
            "t": float(t), "two_sided_p": float(2*stats.t.sf(abs(t), df))}


def rng_for(*tags):
    seed = hashlib.sha256(("nanograv15-single-v1/" + "/".join(map(str, tags))).encode()).digest()
    return np.random.default_rng(int.from_bytes(seed[:8], "big"))


def sensitivity(p, pulsar, name, width, grid, norms, unit, usable, critical):
    """Off-grid injections and a conservative continuous-amplitude upper tolerance bound."""
    rng = rng_for(pulsar, name, width)
    n = INJECTIONS
    summary = {"template": name, "width_d": width, "injections": n,
               "epoch_min_mjd": float(grid.min()), "epoch_max_mjd": float(grid.max())}
    tested = np.flatnonzero(usable)
    if not len(tested):
        return {**summary, "status": "no_testable_template"}
    # Upper confidence bound remains conservative if an actual signal contributes to RSS.
    scale = np.sqrt(p.rss / stats.chi2.ppf(.005, p.dof))
    injections = rng.uniform(grid.min(), grid.max(), n)
    nearest = tested[np.abs(grid[tested, None] - injections).argmin(axis=0)]
    signs = rng.choice([-1., 1.], n)
    bound = np.full(n, np.inf)
    responses = []
    for lo in range(0, n, 32):
        sl = slice(lo, min(lo+32, n))
        # Noise is simulated in coordinates after release-covariance whitening.
        # Any timing/red/DM coefficients cancel under the fixed projection.
        raw = waveform(p, name, injections[sl], width)
        v = p.whiten(raw)
        u = unit[:, nearest[sl]]
        r = remove(p.q, rng.normal(size=v.shape)) * scale
        b = np.sum(u*r, axis=0) * signs[sl]
        c = np.sum(u*v, axis=0)
        r2 = np.maximum(0, np.sum(r*r, axis=0) - b*b)
        v2 = np.maximum(0, np.sum(v*v, axis=0) - c*c)
        rv = np.sum(r*v, axis=0) * signs[sl] - b*c
        k = critical / np.sqrt(p.dof-1)
        margin = c - k*np.sqrt(v2)
        ok = margin > 0
        bnd = np.full(len(c), np.inf)
        bnd[ok] = np.maximum(0, (k*np.sqrt(r2[ok]) - b[ok]) / margin[ok])
        bound[sl] = bnd
        responses.append(np.array([b, c, r2, rv, v2]))
    # For every draw the triangle inequality certifies recovery at all A >= bound.
    # Upper order statistic controls sampling uncertainty of the 96th percentile.
    rank = int(stats.binom.ppf(.995, n, .96)) + 1
    certified = float(np.sort(bound)[rank-1]) if rank <= n else np.inf
    b, c, r2, rv, v2 = np.concatenate(responses, axis=1)
    characteristic = float(np.median(critical * scale / norms[tested]))
    points = np.unique(np.r_[characteristic * np.array([.5, 1, 2, 4]),
                              [certified] if np.isfinite(certified) else []])
    efficiency = []
    for a in points:
        denom = np.sqrt(np.maximum(0, r2 + 2*a*rv + a*a*v2)/(p.dof-1))
        hits = int(np.sum(b + a*c > critical*denom))
        efficiency.append({"amplitude_us": float(a), "recovered": hits})
    return {**summary, "status": "calibrated", "upper_scale": float(scale),
            "power_lower_bound_target": .96, "upper_order_rank": rank,
            "finite_certificate_count": int(np.isfinite(bound).sum()),
            "upper_limit_us_if_no_primary_crossing": certified if np.isfinite(certified) else None,
            "efficiency_nearest_matching_template": efficiency}


def fit_diagnostics(p, unit=None):
    """Compare leverage-adjusted outliers and inter-night structure against projected noise."""
    rng = rng_for("adequacy", len(p.t), p.t.min(), p.t.max())
    leverage = 1 - np.sum(p.q*p.q, axis=1)
    if unit is not None:
        leverage -= unit*unit
    selected = leverage > 1.e-6
    if selected.sum() < 3:
        return {"status": "insufficient_residual_rank", "passes": False}
    _, days = np.unique(np.floor(p.t[selected]), return_inverse=True)
    def metrics(z):
        adjusted = z[selected] / np.sqrt(leverage[selected, None])
        magnitude = np.max(np.abs(adjusted), axis=0) / np.sqrt(np.mean(adjusted**2, axis=0))
        epoch_sums = np.zeros((days.max()+1, z.shape[1]))
        np.add.at(epoch_sums, days, adjusted)
        epoch_sums /= np.sqrt(np.bincount(days)[:, None])
        power = np.maximum(np.sum(epoch_sums**2, axis=0), np.finfo(float).tiny)
        correlation = np.maximum.reduce([np.abs(np.sum(epoch_sums[k:]*epoch_sums[:-k], axis=0))/power
                                         for k in (1, 2, 3)])
        return np.array([magnitude, correlation])
    r = p.ry[:, None]
    if unit is not None:
        r = r - unit[:, None] * (unit @ r)
    observed = metrics(r)
    count = np.zeros(2, dtype=int)
    for _ in range(8):
        noise = remove(p.q, rng.normal(size=(len(p.t), 32)))
        if unit is not None:
            noise -= unit[:, None] * (unit @ noise)
        count += (metrics(noise) >= observed).sum(axis=1)
    prob = (count+1)/257
    return {"status": "tested", "outlier_p": float(prob[0]),
            "lag_structure_p": float(prob[1]), "passes": bool(np.all(prob >= .005))}


def physical_amplitude(name, amplitude):
    if name == "spin_step_or_memory":
        return {"fractional_frequency_step_if_intrinsic": -amplitude * 1.e-6/(YEAR*86400)}
    if name == "spin_slope_step":
        return {"fractional_frequency_derivative_step_per_s": -2*amplitude*1.e-6/(YEAR*86400)**2}
    if name.startswith("dm_"):
        return {"delta_dm_pc_cm3_if_dispersive": amplitude * 1400**2 / 4.148808e9}
    return {}


def row(name, epoch, width, t, amp, se, prob, multiplicity):
    return {"template": name, "epoch_mjd": float(epoch), "width_d": width,
            "amplitude_us": float(amp), "se_us": float(se), "t": float(t),
            "p_local": float(prob), "p_global": min(1., 2*multiplicity*float(prob)),
            **physical_amplitude(name, float(amp))}


def analyze(data, path):
    pulsar = path.stem
    p = fit_pulsar(path)
    grid = epoch_grid(p.t)
    groups = [(s, None) for s in FIXED] + [(s, w) for s in SHAPES for w in WIDTHS]
    trials = len(grid)*len(groups)
    multiplicity = trials*PULSARS
    critical = stats.t.isf(SINGLE_ALPHA / (2*multiplicity), p.dof-1)
    result = {"pulsar": pulsar, "status": "ok", "trials": trials,
              "critical_abs_t": float(critical), "coverage": {
                  "rows": len(p.t), "distinct_days": len(np.unique(np.floor(p.t))),
                  "first_mjd": float(p.t.min()), "last_mjd": float(p.t.max()),
                  "first_search_epoch_mjd": float(grid.min()), "last_search_epoch_mjd": float(grid.max()),
                  "epoch_count": len(grid), "backends": sorted(set(p.backend.tolist())),
                  "dm_restored": p.dm_restored, "dm_retained": p.dm_retained},
              "noise": p.noise, "noise_diagnostics": fit_diagnostics(p), "peaks": [],
              "primary_crossings": [], "survivors": [], "sensitivity": []}
    for name, width in groups:
        raw = waveform(p, name, grid, width)
        norm, unit, t, amp, se = score(p, p.whiten(raw))
        usable = measurable(p, raw, norm)
        prob = 2*stats.t.sf(abs(t), p.dof-1)
        indices = np.flatnonzero(usable)
        if len(indices):
            best = indices[np.argmax(abs(t[indices]))]
            result["peaks"].append(row(name, grid[best], width, t[best], amp[best], se[best],
                                       prob[best], multiplicity))
        crossing = indices[abs(t[indices]) >= critical]
        for i in crossing:
            result["primary_crossings"].append(row(name, grid[i], width, t[i], amp[i], se[i],
                                                   prob[i], multiplicity))
        calibration = sensitivity(p, pulsar, name, width, grid, norm, unit, usable, critical)
        calibration["primary_crossing_count"] = len(crossing)
        reasons = []
        if len(crossing):
            reasons.append("primary_crossing_in_family_and_width")
        if not result["noise_diagnostics"]["passes"]:
            reasons.append("noise_model_diagnostic")
        if name.startswith("dm_") and not p.dm_retained:
            reasons.append("dm_retention_unverified")
        if calibration.get("upper_limit_us_if_no_primary_crossing") is None:
            reasons.append("no_finite_upper_certificate")
        calibration["exclusion_withheld"] = reasons
        calibration["reported_upper_limit_us"] = (None if reasons else
                                                  calibration["upper_limit_us_if_no_primary_crossing"])
        result["sensitivity"].append(calibration)

    checks_data = {}
    for dataset in ("wb", "epta"):
        fn = data / dataset / path.name
        if not fn.is_file():
            checks_data[dataset] = "not_available"
            continue
        try:
            checks_data[dataset] = fit_pulsar(fn, fit_soft=False)
        except Exception as exc:
            checks_data[dataset] = "reduction_unusable: " + str(exc)
    result["available_checks"] = {key: "loaded" if isinstance(val, Projection) else val
                                  for key, val in checks_data.items()}
    for event in result["primary_crossings"]:
        name = event["template"]
        raw = waveform(p, name, event["epoch_mjd"], event["width_d"]).ravel()
        soft_template = p.soft_whiten(raw)
        ss = float(soft_template @ soft_template)
        soft_estimate = float(soft_template @ p.soft_whiten(p.y) / ss) if ss > 0 else 0.
        spectra = chromatic_fit(p, raw, name)
        resolved = (spectra.get("expected_exponent_partial_p", 1.) < .01
                    and spectra.get("expected_exponent_partial_t", 0.) * event["amplitude_us"] > 0)
        event["spectrum"] = "conditional_component_supported" if resolved else "ambiguous"
        checks = {"red_dm_reml_amplitude_us": soft_estimate, "leave_one_day": delete_visit(p, raw),
                  "joint_frequency_dependence": spectra}
        for dataset, other in checks_data.items():
            checks[dataset] = fixed_check(other, event) if isinstance(other, Projection) else {
                "status": other}
        reasons = []
        if soft_estimate * event["amplitude_us"] <= 0 or not p.noise.get("reml_converged", False):
            reasons.append("red_dm_reml_check")
        leave = checks["leave_one_day"]
        if (leave["status"] != "tested" or leave.get("two_sided_p", 1.) * multiplicity > SINGLE_ALPHA
                or leave.get("amplitude_us", 0.) * event["amplitude_us"] <= 0):
            reasons.append("leave_one_day_check")
        if checks["wb"].get("directional_p", 1.) > .05:
            reasons.append("wideband_same_observation_check")
        if checks["epta"].get("directional_p", 1.) > .05/max(1, len(result["primary_crossings"])):
            reasons.append("independent_archive_check")
        if name.startswith("dm_"):
            if not all([p.dm_retained, checks["wb"].get("dm_retained", False),
                        checks["epta"].get("dm_retained", False)]):
                reasons.append("dm_retention_unverified")
            wb = checks_data["wb"]
            dm_check = wideband_dm_check(wb, event) if isinstance(wb, Projection) else {"status": wb}
            checks["wideband_direct_dm"] = dm_check
            if dm_check.get("directional_p", 1.) > .05:
                reasons.append("wideband_direct_dm_check")
        if not reasons:
            z = p.whiten(raw)
            checks["post_event_noise_diagnostics"] = fit_diagnostics(p, z / np.linalg.norm(z))
            if not checks["post_event_noise_diagnostics"]["passes"]:
                reasons.append("post_event_noise_diagnostic")
        event["checks"] = checks
        event["unpassed"] = reasons
        if not reasons:
            result["survivors"].append(event)
    result["peak_neighborhoods"] = group_crossings(result["primary_crossings"])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    args = parser.parse_args()
    if args.shards < 1 or not 0 <= args.shard < args.shards:
        parser.error("--shard must be within [0,--shards)")
    out = args.out / "single"
    out.mkdir(parents=True, exist_ok=True)
    paths = sorted((args.data/"nb").glob("*.npz"))
    if len(paths) > PULSARS:
        raise ValueError("more primary targets than prespecified family allocation")
    for path in paths[args.shard::args.shards]:
        try:
            result = analyze(args.data, path)
        except Exception as exc:
            result = {"pulsar": path.stem, "status": "error", "error": str(exc),
                      "traceback": traceback.format_exc()}
        target = out / (path.stem + ".json")
        temp = target.with_suffix(".tmp")
        temp.write_text(json.dumps(result, allow_nan=False, separators=(",", ":"))+"\n")
        temp.replace(target)
        print(path.stem, result["status"], flush=True)


if __name__ == "__main__":
    main()
