"""Inject widely separated lines into real calibrated spectra at fixed stellar templates."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np

from laser_search import inject_line, scan_target
from recovery import _matches

MULTIPLIERS = (8., 16., 32., 64., 128.)


def local_response_fwhm(e, bary):
    # Mix individually flux-normalized columns as in the detector.
    w = e.wave
    j = np.searchsorted(w, bary, side="right") - 1
    m = e.resolution.shape[0]
    h = m // 2
    if j-h < 0 or j+h+1 >= len(w):
        return None
    phase = (bary-w[j]) / (w[j+1]-w[j])
    widths = np.gradient(w)[j-h:j+h+2]
    p = np.zeros(m+1)
    for i, frac in ((0, 1-phase), (1, phase)):
        col = np.asarray(e.resolution[:, j+i], float)
        area = col @ widths[i:i+m]
        if not np.isfinite(area) or area <= 0:
            return None
        p[i:i+m] += frac * col / area
    mass = np.maximum(0, p)*widths
    wc = w[j-h:j+h+2]
    mean = np.average(wc, weights=mass)
    return float(2.35482*np.sqrt(np.average((wc-mean)**2, weights=mass)))


def noise_scale(e):
    """Median per-column uncertainty of optimal unit-area flux, before nuisance fitting."""
    w = np.asarray(e.wave)
    r = np.asarray(e.resolution, float)
    m,n = r.shape
    h = m//2
    j = np.arange(h, n-h)
    rows = j[:, None]+np.arange(-h,h+1)
    ww = np.asarray(e.ivar, float)[rows]
    col = r[:,j].T
    use = (np.asarray(e.mask)[rows] == 0).all(1) & (ww > 0).all(1)
    area = (col*np.gradient(w)[rows]).sum(1)
    norm = (ww*col**2).sum(1)
    use &= np.isfinite(area) & (area > 0) & np.isfinite(norm) & (norm > 0)
    if not np.any(use):
        return None
    return float(np.median(area[use] / np.sqrt(norm[use])))


def separated_groups(indices, rest_waves, exposures):
    # Leave a gap large enough for disjoint injection support, baseline windows,
    # and peak/centroid decisions. Inputs are actual (not artificially concatenated) grids.
    separation = max((3*e.resolution.shape[0] + 20)*np.max(np.diff(e.wave))/e.doppler
                     for e in exposures)
    groups = []
    for i in sorted(indices, key=lambda k: rest_waves[k]):
        for group in groups:
            if rest_waves[i] - rest_waves[group[-1]] > separation:
                group.append(i)
                break
        else:
            groups.append([i])
    return groups


def experiment(exposures, expected, n_trials, seed):
    exposures = tuple(exposures)
    base = scan_target(exposures, expected)
    lo = max(e.wave[0]/e.doppler for e in exposures)
    hi = min(e.wave[-1]/e.doppler for e in exposures)
    rest = np.random.default_rng(seed).uniform(lo, hi, n_trials)
    scales = [noise_scale(e) for e in exposures]
    if any(s is None for s in scales):
        raise ValueError("no valid resolution columns for sensitivity")
    fluxes = [float(m * max(scales)) for m in MULTIPLIERS]
    positions = []
    selected = []
    for i, w in enumerate(rest):
        eligible = any(a <= w <= b for a,b in base.coverage)
        confused = any(_matches(candidate,w) for candidate in base.candidates)
        reason = "outside_joint_coverage" if not eligible else ("baseline_coincidence" if confused else None)
        fwhm = [local_response_fwhm(e, w*e.doppler) for e in exposures] if eligible else []
        positions.append(dict(rest_A=float(w), eligible=eligible, fwhm_A=fwhm,
                              reasons=[reason]*len(fluxes), fitted=[None]*len(fluxes)))
        if reason is None:
            selected.append(i)
    groups = separated_groups(selected, rest, exposures)
    for level, flux in enumerate(fluxes):
        for group in groups:
            injected = list(exposures)
            for i in group:
                injected = [inject_line(e, rest[i], flux) for e in injected]
            result = scan_target(injected, expected)
            for i in group:
                near = [c for c in result.candidates if _matches(c, rest[i])]
                if near:
                    c = min(near, key=lambda c: abs(c.rest_A-rest[i]))
                    by_id = {h.expid:h for h in c.hits}
                    positions[i]["fitted"][level] = [[float(by_id[e.expid].observed_A),
                                                       float(by_id[e.expid].rest_fwhm_A*e.doppler)]
                                                      for e in exposures]
                    positions[i]["reasons"][level] = "recovered"
                else:
                    positions[i]["reasons"][level] = "missed_detection"
    return dict(seed=int(seed), range_rest_A=[float(lo),float(hi)], flux17=fluxes,
                resolution_noise17=scales, trials=positions,
                epoch=[dict(expid=int(e.expid), doppler=float(e.doppler), rv_error=float(e.velocity_error_kms))
                       for e in exposures])


def run(cache, survey, program, hp, audit_file, output, n_trials, partition=0, partitions=1):
    from desi_shard import targets
    with open(audit_file) as f:
        audit = {int(r["targetid"]): r for r in csv.DictReader(f)
                 if r["all_assigned"]=="True" and r["star_class_consistent"]=="True"
                 and r["ambiguous_gaia"]=="False"}
    if not 0 <= partition < partitions:
        raise ValueError("invalid partition index or number of partitions")
    suffix = f"-part{partition}of{partitions}" if partitions > 1 else ""
    dest = Path(output) / f"injections-{survey}-{program}-{hp}{suffix}.jsonl.gz"
    tmp = dest.with_suffix(".partial")
    with gzip.open(tmp, "wt", compresslevel=5) as file:
        for meta, arms in targets(cache, survey, program, hp):
            tid = meta["targetid"]
            if tid not in audit or tid % partitions != partition:
                continue
            expected = tuple(int(v) for v in audit[tid]["expected_expids"].split("|"))
            for ai,(arm, records) in enumerate(sorted(arms.items())):
                seed = ((tid << 2) + ai)
                res = experiment(records, expected, n_trials, seed)
                print(json.dumps(dict(targetid=tid, arm=arm, **res), allow_nan=False, separators=(',', ':')),
                      file=file, flush=True)
            print("complete", tid, flush=True)
    tmp.replace(dest)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache")
    p.add_argument("survey")
    p.add_argument("program")
    p.add_argument("healpix", type=int)
    p.add_argument("audit")
    p.add_argument("output")
    p.add_argument("--trials", type=int, default=160)
    p.add_argument("--partition", type=int, default=0)
    p.add_argument("--partitions", type=int, default=1)
    a = p.parse_args()
    run(a.cache,a.survey,a.program,a.healpix,a.audit,a.output,a.trials,a.partition,a.partitions)
