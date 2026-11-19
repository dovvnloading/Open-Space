#!/usr/bin/env python3
"""Evaluate diagnostics from *real paired* per-source stellar mass/redshift draws.

Input CSV REQUIRED columns:
    source_id,sed_model,z,log10_mstar,weight
Each row must be a JOINT draw from one SED fit (never independently paired
marginal mass and redshift quantiles). sed_model distinguishes independent fits.

Requires halo_bounds.json generated from halo_bound.py, including a full
lambda(>M) curve for each z bin. Do not supply or interpret simulated toy
chains as observations. Eddington bias and survey selection are NOT handled:
this is a conditional screening statistic, not a LCDM p-value.
"""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

FIELDS = ("source_id", "sed_model", "z", "log10_mstar", "weight")


def evaluate(rows, bounds):
    bins = bounds["bins"]
    grid = np.log10(np.asarray(bounds["halo_mass_grid_msun"], float))
    if np.any(np.diff(grid) <= 0):
        raise ValueError("halo mass grid must increase")
    by_target = defaultdict(list)
    for item in rows:
        if any(x not in item for x in FIELDS):
            raise ValueError("Input requires paired JOINT posterior samples: " + ", ".join(FIELDS))
        v = {k: item[k] for k in FIELDS}
        for col in ("z", "log10_mstar", "weight"):
            v[col] = float(v[col])
        if not all(math.isfinite(v[col]) for col in ("z", "log10_mstar", "weight")) or v["weight"] < 0:
            raise ValueError("invalid joint posterior sample")
        by_target[(v["source_id"], v["sed_model"])].append(v)
    out = []
    fb = bounds["f_baryon"]
    for (source_id, model_name), draws in sorted(by_target.items()):
        weights = np.array([x["weight"] for x in draws], float)
        n_eff = float(weights.sum()**2 / (weights @ weights)) if weights.sum() else 0.
        if n_eff <= 0:
            raise ValueError("Nonpositive normalized weights for " + source_id)
        weights /= weights.sum()
        r = {"source_id":source_id,"sed_model":model_name,"n_draws":len(draws),
             "effective_sample_size":n_eff,
             "p_z_gt6":float(sum(w for x, w in zip(draws, weights) if x["z"] > 6)),
             "status":"diagnostic only; NOT Eddington-/selection-/AGN-corrected"}
        for hm in ("sheth99", "tinker08"):
            max95 = 0.
            max_expect_one = 0.
            p_no_capable_halo = 0.
            z_in_covered = 0.
            for row, w in zip(draws, weights):
                z, logm = row["z"], row["log10_mstar"]
                if not (6 < z <= 15):
                    continue
                match = next((b for b in bins if b["z_min_exclusive"] < z <= b["z_max_inclusive"]), None)
                if match is None:
                    continue
                z_in_covered += w
                h = match["models"][hm]
                quant = h["max_stellar_baryon_quantiles_log10"]["p950"]
                if quant is not None and logm > quant:
                    max95 += w
                if logm > h["log10_M_star_baryon_expected_one"]:
                    max_expect_one += w
                # Conditional probability that a Poisson sample of the NOMINAL
                # volume contains zero halos capable of hosting this draw.
                # Not the posterior probability LCDM is wrong.
                lambdas = np.asarray(h["lambda_expected_count_by_halo_mass_grid"], float)
                if np.any(np.diff(lambdas) > 1e-6):
                    raise ValueError("nonmonotone abundance")
                logreq = logm - math.log10(fb)
                expected = float(np.interp(logreq, grid, lambdas,
                                           left=lambdas[0], right=lambdas[-1]))
                p_no_capable_halo += w * math.exp(-max(0., expected))
            r[hm] = {"p_joint_z_gt6_above_p95_max_halo_baryon":float(max95),
                     "p_joint_z_gt6_above_expected_one_baryon":float(max_expect_one),
                     "posterior_mean_poisson_zero_capable_halo_diagnostic":float(p_no_capable_halo),
                     "posterior_prob_z_in_modeled_bins":float(z_in_covered)}
        out.append(r)
    return {"method":"paired SED posterior draws; no independently permuted marginals",
            "limitations":["No Eddington bias or sample selection likelihood",
                           "Approximate nominal-area Poisson halo maximum",
                           "No mass-definition conversion",
                           "No low-z, dwarf, AGN, MIRI, blending, or spectroscopic vetting in this script",
                           "No discovery significance or LCDM rejection inference"],
            "results":out}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--joint-draws", type=Path, required=True)
    p.add_argument("--halo-bounds", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args=p.parse_args()
    with args.joint_draws.open(newline="") as f:
        rows=list(csv.DictReader(f))
    if not rows:
        raise SystemExit("No real joint posterior draws supplied; no exceedance statistic.")
    bounds=json.loads(args.halo_bounds.read_text())
    response=evaluate(rows, bounds)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(response,indent=2)+"\n")
    print(f"wrote {len(response['results'])} source/model posterior diagnostics to {args.out}")


if __name__ == "__main__":
    main()
