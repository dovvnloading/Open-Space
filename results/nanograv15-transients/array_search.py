"""Coherent transient search. Units at the scoring interface are days and microseconds."""
import argparse
import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import beta, chi2, norm

from search import fit_pulsar


YEAR = 365.25
ALPHA = 0.005
SCALE_BUDGET = 0.0005
STEP = 30
SKIES = 48
NF = 60
COMMON_A = 6e-15
MIN_PSR = 8
SHAPES = ("step", "ramp", "gaussian30", "gaussian120")
OFFSETS = (-12., -6., 0., 6., 12.)
NINJ = 256
SEED = 150681


def waveform(shape, dt):
    if shape == "step":
        return (dt >= 0).astype(float)
    if shape == "ramp":
        return np.maximum(dt, 0) / YEAR
    width = float(shape.removeprefix("gaussian"))
    return np.exp(-0.5 * (dt / width) ** 2)


def vector(ra, dec):
    return np.array([np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)])


def sky_grid():
    k = np.arange(SKIES) + 0.5
    ra = 2*np.pi * np.mod(k*(np.sqrt(5)-1)/2, 1)
    dec = np.arcsin(1 - 2*k/SKIES)
    return np.array([vector(a, d) for a, d in zip(ra, dec)])


def antenna(xyz, source):
    # source points from Earth to the sky position, opposite the propagation vector.
    ra = np.arctan2(source[1], source[0])
    east = np.array([-np.sin(ra), np.cos(ra), 0.])
    north = np.cross(source, east)
    a, b = xyz @ east, xyz @ north
    denom = 1 - np.clip(xyz @ source, -1, 1)
    return np.column_stack([0.5*(a*a-b*b), a*b]) / np.maximum(denom[:, None], 1e-14)


def common_angular(xyz):
    dot = np.clip(xyz @ xyz.T, -1, 1)
    x = (1 - dot) / 2
    hd = 0.5 - 0.25*x + 1.5*x*np.log(np.maximum(x, 1e-300))
    hd += 0.5 * np.eye(len(dot))
    return hd + np.ones_like(hd) + dot


def common_basis(t, start, span):
    harmonics = np.arange(1, NF+1)
    angle = 2*np.pi*(t[:, None]-start)*harmonics/span
    # P_r(f) df in us^2, characteristic strain at 1/Julian year.
    sigma = COMMON_A * (YEAR*86400)*1e6 / np.sqrt(12*np.pi**2)
    sigma *= np.sqrt(YEAR/span) * (harmonics*YEAR/span)**(-13/6)
    return np.column_stack([sigma*np.sin(angle), sigma*np.cos(angle)])


def grid(files):
    bounds = []
    for path in files:
        with np.load(path, allow_pickle=False) as f:
            t = f["t"]
            bounds.append((np.min(t), np.max(t)))
    start, stop = min(a for a, _ in bounds), max(b for _, b in bounds)
    epochs = np.arange(STEP*np.ceil(start/STEP), stop, STEP)
    return epochs, float(start), max(float(stop-start), 1.)


def support(t, epochs, shape):
    days = np.unique(np.floor(t))
    ii = np.searchsorted(days, epochs)
    near = 3*float(shape[8:]) if shape.startswith("gaussian") else 180.
    left = np.searchsorted(days, epochs-near)
    right = np.searchsorted(days, epochs+near)
    return (ii >= 6) & (len(days)-ii >= 6) & (ii>left) & (right>ii)


@dataclass
class Pulsar:
    name: str
    xyz: np.ndarray
    dof: int
    rss: float
    z: np.ndarray
    weight: np.ndarray
    common: np.ndarray
    inj_cross: np.ndarray
    inj_rss: np.ndarray
    inj_data: np.ndarray


def reduce_one(path, epochs, start, span, injections):
    p = fit_pulsar(path, fit_soft=False)
    if p.dof < 12:
        raise ValueError(f"only {p.dof} residual degrees of freedom")
    y = p.whiten(p.y)
    red = p.whiten(common_basis(p.t, start, span))
    dims = (len(SHAPES), len(epochs))
    weight, z = np.zeros(dims), np.zeros(dims)
    common = np.zeros((*dims, 2*NF))
    inj_cross, inj_rss, inj_data = (np.zeros((*dims, len(OFFSETS))) for _ in range(3))
    for j, shape in enumerate(SHAPES):
        use = np.flatnonzero(support(p.t, epochs, shape))
        for pos in range(0, len(use), 24):
            k = use[pos:pos+24]
            dt = p.t[:, None] - epochs[k]
            if injections:
                # Apply each sparse jitter block only once across all nearby epochs.
                bank = np.column_stack([waveform(shape,dt-shift) for shift in OFFSETS])
                shifted = np.split(p.whiten(bank),len(OFFSETS),axis=1)
                template = shifted[OFFSETS.index(0.)]
            else:
                template = p.whiten(waveform(shape, dt))
            length = np.linalg.norm(template, axis=0)
            ok = np.isfinite(length) & (length > 1e-10)
            k, template, length = k[ok], template[:, ok], length[ok]
            if not len(k):
                continue
            unit = template / length
            weight[j, k] = length
            z[j, k] = y @ unit
            common[j, k, :] = unit.T @ red
            if injections:
                for l, value in enumerate(shifted):
                    injected = value[:,ok]
                    inj_cross[j, k, l] = np.einsum("ij,ij->j", unit, injected)
                    inj_rss[j, k, l] = np.einsum("ij,ij->j", injected, injected)
                    inj_data[j, k, l] = y @ injected
    return Pulsar(path.stem, vector(p.ra,p.dec), int(p.dof), float(y @ y), z, weight,
                  common, inj_cross, inj_rss, inj_data)


def read_records(files, epochs, start, span, injections):
    records, errors = [], []
    for path in files:
        try:
            records.append(reduce_one(path, epochs, start, span, injections))
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            errors.append({"pulsar": path.stem, "reason": str(exc)})
    return records, errors


@dataclass
class Cell:
    names: list
    xyz: np.ndarray
    weight: np.ndarray
    z: np.ndarray
    common: np.ndarray
    factor: np.ndarray
    rss: np.ndarray
    cross: np.ndarray
    energy: np.ndarray
    data: np.ndarray

    def subset(self, keep):
        idx = np.flatnonzero(keep)
        return Cell([self.names[i] for i in idx], self.xyz[idx], self.weight[idx],
                    self.z[idx], self.common[np.ix_(idx, idx)], self.factor[idx],
                    self.rss[idx], self.cross[idx], self.energy[idx], self.data[idx])


def cell(records, j, k, scale_tests, minimum=MIN_PSR):
    active = [p for p in records if p.weight[j,k] > 0]
    if len(active) < minimum:
        return None
    xyz = np.stack([p.xyz for p in active])
    responses = np.stack([p.common[j,k] for p in active])
    covariance = (responses @ responses.T)*common_angular(xyz)
    dof = np.array([p.dof for p in active])
    # Union upper confidence bounds on all white variances. Covariance at the
    # scoring stage includes a separately specified common-process allowance.
    factor = 1/chi2.ppf(SCALE_BUDGET/scale_tests, dof)
    return Cell([p.name for p in active], xyz, np.array([p.weight[j,k] for p in active]),
                np.array([p.z[j,k] for p in active]), covariance, factor,
                np.array([p.rss for p in active]),
                np.stack([p.inj_cross[j,k] for p in active]),
                np.stack([p.inj_rss[j,k] for p in active]),
                np.stack([p.inj_data[j,k] for p in active]))


def patterns(xyz, family, source):
    monopole = np.ones((len(xyz), 1))
    if family == "clock":
        return monopole, xyz
    if family == "ephemeris":
        return xyz, monopole
    return antenna(xyz, source), np.column_stack([monopole, xyz])


@dataclass
class Fit:
    coefficients: np.ndarray
    covariance: np.ndarray
    statistic: float
    p: float


class Contrast:
    """Fixed linear amplitude estimator; uncertainty uses simultaneous scale bounds."""
    def __init__(self, c, family, source, ref=None):
        g, nuisance = (c.weight[:, None]*v for v in patterns(c.xyz, family, source))
        if ref is None:
            ref = cho_factor(c.common + np.eye(len(c.z)), check_finite=False)
        weighted_g = cho_solve(ref, g, check_finite=False)
        weighted_n = cho_solve(ref, nuisance, check_finite=False)
        if np.linalg.matrix_rank(nuisance) != nuisance.shape[1]:
            raise ValueError("singular nuisance sky geometry")
        filtered = weighted_g - weighted_n @ np.linalg.solve(nuisance.T @ weighted_n,
                                                               nuisance.T @ weighted_g)
        response = filtered.T @ g
        if np.linalg.cond(response) > 1e10 or np.linalg.norm(filtered) < 1e-12:
            raise ValueError("unresolved target sky geometry")
        self.matrix = np.linalg.solve(response, filtered.T)
        self.common = self.matrix @ c.common @ self.matrix.T

    def estimate(self, z, variance):
        b = self.matrix @ z
        cov = self.common + (self.matrix * variance) @ self.matrix.T
        stat = max(float(b @ np.linalg.solve(cov, b)), 0.)
        return Fit(b, cov, stat, float(chi2.sf(stat, len(b))))


def ladder(model, z, variance):
    # Evaluate all injected amplitudes together and re-estimate their white scales.
    b = z @ model.matrix.T
    cov = model.common + np.einsum("pi,mi,qi->mpq", model.matrix, variance, model.matrix)
    inverse_b = np.linalg.solve(cov, b[:, :, None])[:, :, 0]
    return b, np.einsum("mp,mp->m", b, inverse_b)


def trial(c, family, source, ref=None):
    return Contrast(c, family, source, ref).estimate(c.z, c.factor*c.rss)


def physical(shape, family):
    if shape != "ramp":
        return 1., "us"
    if family == "memory":
        return 1e-6/(YEAR*86400), "strain"
    return 1., "us/yr"


def amplitude(coeff):
    return float(coeff[0]) if len(coeff) == 1 else float(np.linalg.norm(coeff))


def angle_fields(source):
    if source is None:
        return {"sky_ra_deg": None, "sky_dec_deg": None}
    return {"sky_ra_deg": float(np.degrees(np.arctan2(source[1],source[0])) % 360),
            "sky_dec_deg": float(np.degrees(np.arcsin(source[2])))}


def p_global(p, count):
    return min(1., SCALE_BUDGET + count*p)


def scan(records, epochs, skies):
    rows, cells = [], {}
    for j, shape in enumerate(SHAPES):
        for k, epoch in enumerate(epochs):
            c = cell(records,j,k,len(records))
            if c is None:
                continue
            cells[j,k] = c
            ref = cho_factor(c.common + np.eye(len(c.z)), check_finite=False)
            tests = [("clock", None), ("ephemeris", None)]
            if shape == "ramp":
                tests.extend(("memory", i) for i in range(len(skies)))
            for family, i in tests:
                source = None if i is None else skies[i]
                try:
                    f = trial(c, family, source, ref)
                except (ValueError, np.linalg.LinAlgError):
                    continue
                scale, units = physical(shape, family)
                rows.append(dict(family=family,shape=shape,epoch_mjd=float(epoch),
                                 sky_id=i,**angle_fields(source),n_pulsars=len(c.names),
                                 amplitude=amplitude(f.coefficients)*scale,units=units,
                                 components=(f.coefficients*scale).tolist(),
                                 component_se=(np.sqrt(np.diag(f.covariance))*scale).tolist(),
                                 largest_axis_se=float(np.sqrt(np.linalg.eigvalsh(f.covariance)[-1]))*scale,
                                 statistic=f.statistic,p_local=f.p,j=j,k=k))
    radius95 = {d: np.sqrt(chi2.isf((0.05-SCALE_BUDGET)/max(1,len(rows)),d)) for d in (1,2,3)}
    for row in rows:
        row["p_global"] = p_global(row["p_local"], len(rows))
        row["p_joint"] = min(1.,2*row["p_global"])
        row["upper95_amplitude"] = (abs(row["amplitude"]) +
                                    radius95[len(row["components"])]*row.pop("largest_axis_se"))
    return rows,cells


def clusters(rows):
    # Crowded sky grid entries and adjacent epochs are kept in the trials ledger.
    groups = {}
    for row in rows:
        if row["p_global"] <= ALPHA:
            groups.setdefault((row["family"],row["shape"]), []).append(row)
    peaks = []
    for group in groups.values():
        members = []
        for row in sorted(group, key=lambda r: r["epoch_mjd"]):
            if members and row["epoch_mjd"]-members[-1]["epoch_mjd"] > 1.1*STEP:
                peaks.append(min(members, key=lambda r:r["p_local"]))
                members = []
            members.append(row)
        if members:
            peaks.append(min(members, key=lambda r:r["p_local"]))
    return peaks


def jackknife(c, family, source, z, variances):
    """Return all leave-one-out statistics for each row of (data, variances).

    Deleting one observation from a GLS fit costs a rank-one normal-matrix
    downdate. Invert the reference pulsar covariance once per geometry.
    """
    n = len(c.names)
    target, nuisance = patterns(c.xyz,family,source)
    d = target.shape[1]
    h = c.weight[:,None]*np.column_stack([target,nuisance])
    ref = cho_factor(c.common + np.eye(n),check_finite=False)
    prec = cho_solve(ref,np.eye(n),check_finite=False)
    u = prec @ h
    invdiag = 1 / prec.diagonal()
    fisher = h.T @ u - np.einsum("ia,ib,i->iab",u,u,invdiag)
    rhs = u.T - np.einsum("ia,ik,i->iak",u,prec,invdiag)
    q = np.zeros((n,len(z)))
    good = np.linalg.cond(fisher) < 1e12
    if not good.any():
        return q
    weights = np.linalg.solve(fisher[good],rhs[good])[:, :d, :]
    weights[np.arange(good.sum()),:,np.flatnonzero(good)] = 0.
    coefficients = weights @ z.T
    common = weights @ c.common @ weights.swapaxes(1,2)
    # precompute W_ri W_si once for the entire amplitude ladder.
    square = (weights[:,:,None,:]*weights[:,None,:,:]).reshape(-1,n)
    cov = (square @ variances.T).reshape(len(weights),d,d,len(z)).transpose(0,3,1,2)
    cov = cov + common[:,None,:,:]
    b = coefficients.transpose(0,2,1)
    inverse_b = np.linalg.solve(cov,b[:,:,:,None])[:,:,:,0]
    q[good] = np.einsum("ima,ima->im",b,inverse_b)
    return q


def leave_one(c, family, source):
    try:
        q = jackknife(c,family,source,c.z[None,:],(c.factor*c.rss)[None,:])[:,0]
    except (ValueError,np.linalg.LinAlgError):
        return dict(status="fail", reason="unresolved leave-one-out covariance")
    i = int(np.argmin(q))
    d = patterns(c.xyz,family,source)[0].shape[1]
    worst = float(chi2.sf(q[i],d))
    return dict(status="pass" if worst < 0.05 else "fail",worst_local_p=worst,
                removed_pulsar=c.names[i])


def replication(c, row, source):
    if c is None:
        return dict(status="no overlapping pulsar brackets the event")
    scale,_ = physical(row["shape"],row["family"])
    b = np.array(row["components"])
    u = b/np.linalg.norm(b)
    sky = patterns(c.xyz,row["family"],source)[0] @ u
    local_scale = 1e6*(YEAR*86400) if row["family"] == "memory" else 1.
    delay_units = "us/yr" if row["shape"] == "ramp" else "us"
    variance = c.factor*c.rss
    individual = [dict(pulsar=name, amplitude=float(c.z[i]/c.weight[i]),
                       se=float(np.sqrt(c.common[i,i]+variance[i])/c.weight[i]),
                       expected_amplitude=float(sky[i]*np.linalg.norm(b)*local_scale),
                       units=delay_units) for i,name in enumerate(c.names)]
    report = dict(n_pulsars=len(c.names), pulsars=c.names, per_pulsar=individual)
    if len(c.names) < MIN_PSR:
        report["status"] = "partial coverage; individual fixed-epoch amplitudes only"
        return report
    try:
        f = trial(c, row["family"], source)
    except (ValueError, np.linalg.LinAlgError):
        report["status"] = "unresolved geometry; individual fixed-epoch amplitudes only"
        return report
    # Freeze discovery direction, sign, waveform, epoch and sky.
    a = float(u @ f.coefficients)*scale
    se = float(np.sqrt(u @ f.covariance @ u))*scale
    p = min(1.,float(norm.sf(a/se))+SCALE_BUDGET)
    report.update(status="same-direction support" if p < .05 else "inconclusive",
                  signed_amplitude=a,units=row["units"],se=se,p_fixed_one_sided=p)
    return report


def write_csv(path, rows, fields):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="", encoding="utf-8") as h:
        writer = csv.DictWriter(h, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(v, separators=(",",":")) if isinstance(v,(dict,list)) else v
                             for k,v in row.items()})


def limit_curve(cells, rows, skies):
    rng = np.random.default_rng(SEED)
    count = len(rows)
    curves, sensitivity = [], []
    for j, shape in enumerate(SHAPES):
        families = ["clock","ephemeris"] + (["memory"] if shape == "ramp" else [])
        for family in families:
            valid = [r for r in rows if r["family"] == family and r["shape"] == shape]
            time_rows = {}
            for r in valid:
                time_rows.setdefault(r["k"], []).append(r)
            eligible = sorted(time_rows)
            if not eligible:
                continue
            # Log grid frozen in timing-delay units. For memory convert us/yr to strain.
            amps = np.r_[0.,2.**np.arange(-8.,13.)]
            scale,units = physical(shape,family)
            recover = np.zeros(len(amps), dtype=int)
            primary = np.zeros_like(recover)
            for _ in range(NINJ):
                k = int(rng.choice(eligible))
                c = cells[j,k]
                off = int(rng.integers(len(OFFSETS)))
                source = None
                if family == "memory":
                    true_source = vector(rng.uniform(0,2*np.pi), np.arcsin(rng.uniform(-1,1)))
                    index = [r["sky_id"] for r in time_rows[k]]
                    source = skies[index[np.argmax(skies[index] @ true_source)]]
                    polarization = rng.uniform(0,2*np.pi)
                    pat = antenna(c.xyz,true_source) @ np.array([np.cos(polarization),
                                                                 np.sin(polarization)])
                elif family == "ephemeris":
                    direction = vector(rng.uniform(0,2*np.pi), np.arcsin(rng.uniform(-1,1)))
                    pat = c.xyz @ direction
                else:
                    pat = np.full(len(c.names), rng.choice([-1.,1.]))
                model = Contrast(c,family,source)
                response = model.matrix @ (pat*c.cross[:,off])
                a = amps[:, None] * pat
                z = c.z + a*c.cross[:,off]
                rss = c.rss + 2*a*c.data[:,off] + a*a*c.energy[:,off]
                variances = c.factor*np.maximum(rss, 0.)
                estimates, q = ladder(model,z,variances)
                hits = ((q >= chi2.isf((ALPHA-SCALE_BUDGET)/count, len(response)))
                        & ((estimates @ response) > 0))
                primary += hits
                # Apply the same all-pulsar leave-one-out check as for the real data.
                local_threshold = chi2.isf(0.05, len(response))
                if hits.any():
                    try:
                        subq = jackknife(c,family,source,z[hits],variances[hits])
                        hits[hits] = subq.min(axis=0) > local_threshold
                    except (ValueError, np.linalg.LinAlgError):
                        hits[:] = False
                recover += hits
            # Simultaneous one-sided exact bounds across the amplitude ladder and
            # all nine family/shape combinations (fixed residual realization).
            tail = 0.05/(len(amps)*9)
            selected = None
            for a, n, n_primary in zip(amps, recover, primary):
                lower = float(beta.ppf(tail,n,NINJ-n+1)) if n else 0.
                r = dict(family=family,shape=shape,amplitude=float(a*scale),units=units,
                         n_trials=NINJ,n_recovered=int(n),n_passing_primary=int(n_primary),
                         recovery_fraction=n/NINJ,
                         lower_95=lower)
                curves.append(r)
                if lower >= 0.9 and selected is None:
                    selected = r["amplitude"]
            worst = max(valid, key=lambda r:r["upper95_amplitude"])
            sensitivity.append(dict(family=family,shape=shape,units=units,
                                    first_tested_amplitude_90pct=selected,
                                    worst_eligible_grid_upper95=worst["upper95_amplitude"],
                                    worst_epoch_mjd=worst["epoch_mjd"],worst_sky_id=worst["sky_id"]))
    return curves,sensitivity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    files = sorted((args.input/"nb").glob("*.npz"))
    summary = dict(global_alpha=ALPHA, scale_budget=SCALE_BUDGET, n_pulsars=0,
                   n_trials=0, survivors=[], sensitivity=[], input_failures={},
                   p_global_scope="all eligible array trials",
                   p_joint_scope="both individual and coherent searches")
    fields = ("family shape epoch_mjd sky_id sky_ra_deg sky_dec_deg n_pulsars amplitude units "
              "components component_se statistic p_local p_global p_joint upper95_amplitude").split()
    checks, curves, rows = [], [], []
    if files:
        epochs,start,span = grid(files)
        records,errors = read_records(files,epochs,start,span,True)
        summary.update(n_pulsars=len(records), participating_pulsars=[p.name for p in records])
        summary["input_failures"]["nb"] = errors
        skies = sky_grid()
        rows,cells = scan(records,epochs,skies)
        summary["n_trials"] = len(rows)
        summary["strongest_trial"] = ({k:v for k,v in min(rows,key=lambda r:r["p_local"]).items()
                                       if k in fields} if rows else None)
        peaks = clusters(rows)
        checks_by_dataset = {}
        for dataset in ("wb","epta"):
            paths = [path for path in sorted((args.input/dataset).glob("*.npz"))
                     if path.stem in summary["participating_pulsars"]]
            if not peaks or not paths:
                checks_by_dataset[dataset] = {}
                continue
            # Only the discovery epochs are queried in the validation reductions.
            want = np.unique([r["epoch_mjd"] for r in peaks])
            validation,errors = read_records(paths,want,start,span,False)
            summary["input_failures"][dataset] = errors
            checks_by_dataset[dataset] = {}
            for r in peaks:
                # Restrict geometric validation to the same participating pulsars.
                nb_names = set(cells[r["j"],r["k"]].names)
                subset = [p for p in validation if p.name in nb_names]
                checks_by_dataset[dataset][r["j"],r["k"]] = cell(
                    subset,r["j"],int(np.searchsorted(want,r["epoch_mjd"])),len(paths),minimum=1)
        for row in peaks:
            j,k = row["j"],row["k"]
            c = cells[j,k]
            source = None if row["sky_id"] is None else skies[row["sky_id"]]
            event = {key:row[key] for key in fields}
            event.update(pulsars=c.names, internal_check=leave_one(c,row["family"],source))
            for dataset,label in (("wb","wideband_check"), ("epta","epta_check")):
                event[label] = replication(checks_by_dataset[dataset].get((j,k)),row,source)
            checks.append(event)
            if event["internal_check"]["status"] == "pass":
                summary["survivors"].append(event)
        curves,summary["sensitivity"] = limit_curve(cells,rows,skies)
    summary["status"] = "analyzed" if rows else "inconclusive: insufficient eligible sky geometry"
    summary["threshold_clusters"] = len(checks)
    summary["calibration"] = ("Finite Fourier intrinsic/DM nuisance; additive independent Gaussian "
        "white scales bounded simultaneously with 0.0005 failure probability; 60 common Fourier "
        "harmonics, gamma=13/3, max equivalent strain 6e-15 each for HD, monopole and dipole. "
        "Fixed geometry contrasts with chi-square union bound. Valid conditional on this noise model.")
    summary["upper95_amplitude"] = ("Simultaneous one-sided upper bound on the Euclidean norm of the "
        "response amplitude on the eligible grid, conditional on the fixed waveform and Gaussian "
        "noise model. Worst eligible-grid bounds by family/shape accompany injection completeness.")
    summary["recovery"] = ("256 off-grid injections per family/shape into observed NG narrowband "
        "residuals, uniform eligible 30-day grid epochs with offsets -12,-6,0,6,12 d and isotropic "
        "directions/polarizations; nearest-cell family-wise crossing with aligned amplitude and "
        "all leave-one-out tests passing. Recompute variance bounds after each injection. "
        "Completeness is conditional on this residual realization and precedes external characterization.")
    write_csv(args.output/"array_trials.csv.gz",rows,fields)
    write_csv(args.output/"array_checks.csv",checks,fields+["pulsars","internal_check",
                                                            "wideband_check","epta_check"])
    write_csv(args.output/"array_limits.csv",curves,
              "family shape amplitude units n_trials n_recovered n_passing_primary recovery_fraction lower_95".split())
    (args.output/"array_summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
