"""Post-selection fixed-template checks of narrowband array fits against local dispersion."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

from array_search import (MIN_PSR, NF, OFFSETS, SHAPES, Pulsar, cell, common_basis,
                          grid, replication, sky_grid, support, vector, waveform, write_csv)
from search import Projection, basis, remove

LEVERAGE = .01
VARIANTS = ("dm_free_all", "dm_free_arecibo", "dm_free_gbt",
            "matched_raw_all", "matched_raw_arecibo", "matched_raw_gbt")


def physical_site(backend):
    name = str(backend).split("/", 1)[0].lower()
    if name in ("ao", "arecibo") or name.startswith("arecibo_"):
        return "arecibo"
    if name in ("gbt", "greenbank", "green_bank"):
        return "gbt"
    return name


def intercept_weights(err, chrom, group, jitter):
    """GLS weights for the coefficient of 1 in the design [1, (1400/f)^2]."""
    if len(err) < 2:
        return None
    x = np.column_stack((np.ones(len(err)), chrom))
    inverse_variance = 1 / err**2
    inverse_x = inverse_variance[:, None] * x
    # Apply (diag(err²) + sum u_g u_g^T)^-1. ECORR support sets are disjoint.
    for g in np.unique(group[group >= 0]):
        ix = np.flatnonzero(group == g)
        u = jitter[ix]
        v = inverse_variance[ix] * u
        inverse_x[ix] -= np.outer(v, v @ x[ix]) / (1 + u @ v)
    gram = x.T @ inverse_x
    leverage = 1 - gram[0,1]**2 / (gram[0,0]*gram[1,1])
    if not np.isfinite(leverage) or leverage < LEVERAGE:
        return None
    first = np.linalg.solve(gram, [1., 0.])
    w = inverse_x @ first
    if not np.allclose(w @ x, [1., 0.], rtol=0, atol=1e-9) or first[0] <= 0:
        raise ValueError("failed intercept or DM cancellation identity")
    return w, np.sqrt(first[0]), leverage


def visits(data, sites):
    t, err = data["t"], data["err"]
    chrom = (1400 / data["freq"])**2
    group, jitter = data["ecorr_group"], data["ecorr_us"]
    day = np.floor(t).astype(int)
    # Keep an ECORR block entirely within one GLS fit, including near midnight.
    for group_id in np.unique(group[group >= 0]):
        ix = np.flatnonzero(group == group_id)
        if len(np.unique(sites[ix])) != 1:
            raise ValueError("ECORR group straddles observatories")
        day[ix] = int(np.floor(np.median(t[ix])))
    members = {}
    for i, key in enumerate(zip(sites, day)):
        members.setdefault(key, []).append(i)
    row_index, col_index, weights = [], [], []
    times, errors, site_keys, leverage = [], [], [], []
    census = {}
    for (site, _), indices in members.items():
        count = census.setdefault(site, {"observed_site_days": 0, "usable": 0})
        count["observed_site_days"] += 1
        ix = np.asarray(indices)
        fit = intercept_weights(err[ix], chrom[ix], group[ix], jitter[ix])
        if fit is None:
            continue
        w, se, lever = fit
        row_index.extend([len(times)] * len(ix))
        col_index.extend(indices)
        weights.extend(w)
        times.append(np.average(t[ix], weights=err[ix]**-2))
        errors.append(se)
        site_keys.append(site)
        leverage.append(lever)
        count["usable"] += 1
    L = csr_matrix((weights, (row_index, col_index)), shape=(len(times), len(t)))
    return L, np.asarray(times), np.asarray(errors), np.asarray(site_keys), census, leverage


def compressed_projection(data, matrix, t, err):
    """Propagate the published timing Jacobian and system phase offsets through the intercept."""
    if len(t) < 16 or np.ptp(t) <= 0:
        raise ValueError("fewer than 16 useful site-days or zero span")
    rt = data["t"]
    x = (rt - np.mean(t)) / np.ptp(t)
    be = data["backend"]
    raw_design = np.column_stack((data["design"], be[:,None] == np.unique(be), x, x*x))

    def propagate(columns):
        compressed = matrix @ columns
        # Signed intercept weights exactly cancel a pure dispersion column.
        # Guard round-off before normalizing columns for rank determination.
        scale = np.max(abs(matrix) @ np.abs(columns), axis=0)
        keep = np.max(np.abs(compressed), axis=0) > 1e-11 * scale
        return compressed[:,keep] / err[:,None]

    design = propagate(raw_design)
    qtime = basis(design)
    # A free DM per site-day has already been profiled out. Achromatic Fourier
    # amplitudes remain free. Use the discovery rule for harmonic spacing/cap.
    budget = len(t) - qtime.shape[1] - 14
    modes = min(20, max(1, int(np.ptp(t) // 180)), max(0, budget // 4))
    phase = 2*np.pi*x[:,None] * np.arange(1, modes+1)
    if modes:
        design = np.column_stack((design, propagate(np.column_stack((np.sin(phase), np.cos(phase))))))
    q = basis(design)
    def project(raw):
        return remove(q, (matrix @ raw) / (err if raw.ndim == 1 else err[:,None]))
    return project, len(t)-q.shape[1]


def make_record(name, ra, dec, t, raw_t, raw_y, project, dof, epochs, start, span):
    if dof < 12:
        raise ValueError(f"only {dof} residual degrees of freedom")
    y = project(raw_y)
    red = project(common_basis(raw_t, start, span))
    dims = (len(SHAPES), len(epochs))
    weight, z = np.zeros(dims), np.zeros(dims)
    common = np.zeros((*dims, 2*NF))
    for j, shape in enumerate(SHAPES):
        use = np.flatnonzero(support(t, epochs, shape))
        for pos in range(0, len(use), 24):
            k = use[pos:pos+24]
            template = project(waveform(shape, raw_t[:, None] - epochs[k]))
            norms = np.linalg.norm(template, axis=0)
            keep = np.isfinite(norms) & (norms > 1e-10)
            idx = k[keep]
            unit = template[:,keep] / norms[keep]
            weight[j,idx] = norms[keep]
            z[j,idx] = y @ unit
            common[j,idx] = unit.T @ red
    zeros = np.zeros((*dims, len(OFFSETS)))
    return Pulsar(name, vector(ra,dec), dof, float(y @ y), z, weight, common, zeros, zeros, zeros)


def read_controls(path, epochs, start, span):
    with np.load(path, allow_pickle=False) as src:
        data = {key: src[key] for key in ("t", "y", "err", "freq", "design", "backend",
                                           "ecorr_group", "ecorr_us", "ra", "dec")}
    sites = np.array([physical_site(be) for be in data["backend"]])
    L, t, err, obs, census, leverage = visits(data, sites)
    report = {"pulsar":path.stem, "raw_rows":len(sites), "site_days":census,
              "minimum_accepted_leverage": min(leverage) if leverage else None,
              "controls":{}}
    records = {}
    for variant in VARIANTS:
        site = variant.rsplit("_",1)[-1]
        try:
            keep = np.ones(len(t),dtype=bool) if site == "all" else obs == site
            if variant.startswith("matched_raw"):
                # Same sites, days and arrival-time measurements as the
                # DM-free analysis. Remaining differences expose the cost of
                # fitting a separate DM at each site-day and compressing it.
                raw_rows = np.unique(L[keep].indices)
                cut = {key: val[raw_rows] if val.ndim else val for key,val in data.items()}
                p = Projection(cut, fit_soft=False)
                used, raw_t, y, project, dof = p.t, p.t, p.y, p.whiten, p.dof
            else:
                used = t[keep]
                project, dof = compressed_projection(data, L[keep], used, err[keep])
                raw_t, y = data["t"], data["y"]
            rec = make_record(path.stem, float(data["ra"]), float(data["dec"]), used, raw_t,
                              y, project, dof, epochs, start, span)
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            report["controls"][variant] = {"status":"unavailable", "reason":str(exc)}
            continue
        records[variant] = rec
        report["controls"][variant] = {"status":"available", "n_rows":len(used),
                                      "mjd_min":float(min(used)), "mjd_max":float(max(used)),
                                      "dof":int(dof)}
    return records, report


def calibration_check():
    # A time-independent reference delay and DM column recover exactly, with
    # variance equal to the full dense block-covariance calculation.
    freq = np.array([780., 830., 1400., 1510.])
    chrom = (1400/freq)**2
    err = np.array([.75, .65, .8, .9])
    ids = np.array([1,1,2,2])
    jitter = np.array([.4,.4,.3,.3])
    w, se, _ = intercept_weights(err, chrom, ids, jitter)
    C = np.diag(err**2) + (ids[:,None] == ids) * np.outer(jitter,jitter)
    np.testing.assert_allclose(w @ (3.4-112*chrom), 3.4, atol=1e-9)
    np.testing.assert_allclose(w @ C @ w, se**2, rtol=1e-10)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", required=True,type=Path)
    p.add_argument("--evidence", required=True,type=Path)
    p.add_argument("--out", required=True,type=Path)
    args = p.parse_args()
    calibration_check()
    args.out.mkdir(parents=True, exist_ok=True)
    summary = json.loads((args.evidence/"array/array_summary.json").read_text())
    events = summary["survivors"]
    want = np.unique([event["epoch_mjd"] for event in events])
    paths = sorted((args.data/"nb").glob("*.npz"))
    records = {v:{} for v in VARIANTS}
    ledger = {}
    if events:
        _, start, span = grid(paths)
        for path in paths:
            try:
                variants, ledger[path.stem] = read_controls(path, want, start, span)
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
                ledger[path.stem] = {"status":"unavailable","reason":str(exc)}
                continue
            for variant, record in variants.items():
                records[variant][path.stem] = record
    rows, full = [], []
    for ordinal, event in enumerate(events, 1):
        ident = {k:event[k] for k in ("family","shape","epoch_mjd","sky_id","units")}
        ident["survivor_ordinal"] = ordinal
        source = None if ident["sky_id"] is None else sky_grid()[ident["sky_id"]]
        k, j = int(np.searchsorted(want, event["epoch_mjd"])), SHAPES.index(event["shape"])
        for variant, catalog in records.items():
            site = variant.rsplit("_",1)[-1]
            companion = records[("matched_raw_" if variant.startswith("dm_free_") else "dm_free_") + site]
            base = [n for n in event["pulsars"] if n in catalog and catalog[n].weight[j,k] > 0]
            pool = [catalog[n] for n in base if n in companion and companion[n].weight[j,k] > 0]
            c = cell(pool,j,k,max(len(paths),1),minimum=1)
            checked = replication(c, event, source)
            detail = {**ident, "control":variant, "discovery_pulsars":len(event["pulsars"]),
                      "eligible_before_membership_match":len(base),
                      "n_pulsars":0, **checked}
            reference = float(np.linalg.norm(event["components"]))
            detail["discovery_amplitude_along_fixed_axis"] = reference
            if detail.get("se", 0.) > 0 and reference > 0:
                detail["expected_snr_at_discovery_amplitude"] = reference / detail["se"]
                detail["signed_amplitude_ratio"] = detail["signed_amplitude"] / reference
            full.append(detail)
            rows.append({k:v for k,v in detail.items() if k != "per_pulsar"})
    fields = ("survivor_ordinal family shape epoch_mjd sky_id units control discovery_pulsars "
              "eligible_before_membership_match n_pulsars "
              "status signed_amplitude se p_fixed_one_sided discovery_amplitude_along_fixed_axis "
              "expected_snr_at_discovery_amplitude signed_amplitude_ratio pulsars").split()
    write_csv(args.out/"fixed_controls.csv", rows, fields)
    notes = {"description":"Post-selection, same-data diagnostics evaluated at the narrowband discovery "
             "epoch, shape and signed response. Local one-sided p values are descriptive.",
             "leverage_rule":"At each site-day require 1 - rho^2 >= 0.01 for GLS columns 1 and (1400 MHz/f)^2; "
             "maximum permitted intercept variance inflation is 100. ECORR groups stay intact.",
             "comparison":"Both members of each raw/intercept pair use the identical participating set at each "
             "event, on identical site-days. Their residual degrees of freedom, red harmonic budgets and "
             "noise precision can differ. Amplitude changes alone do not identify a dispersive origin.",
             "minimum_coherent_pulsars":MIN_PSR, "events":len(events), "control_fits":len(full),
             "per_pulsar":ledger}
    (args.out/"coverage.json").write_text(json.dumps(notes,indent=2,allow_nan=False)+"\n")
    (args.out/"fixed_controls.json").write_text(json.dumps(full,allow_nan=False,separators=(",",":"))+"\n")


if __name__ == "__main__":
    main()
