"""Microlensing-timescale sensitivity: likelihood competition with independent variability.

For each bin width and *fixed* slow-microlensing timescale, refit remaining GP
hyperparameters at delta=210, then scan delay at the fitted hyperparameters.
An independent-images covariance with its own fitted nuisance is a competing
model, not a p-value or a Bayes factor. All values are conditional on choices.
"""
import argparse
import json
import numpy as np
from scipy.optimize import minimize
from gaussian_process import DelayedGP, observations, optimize, INITIAL, BOUNDS


def stress(width, ell):
    gp = DelayedGP(observations("ztf", width=width))
    init = INITIAL.copy()
    init[4] = np.log(ell)
    bounds = BOUNDS.copy()
    bounds[4] = np.log([ell, ell])
    opt = minimize(lambda x:gp.score(210, x),init,method="L-BFGS-B",bounds=bounds,
                   options={"maxiter":60,"ftol":2e-6})
    grid = range(-100,601,5)
    profiles = {d:gp.score(d, opt.x) for d in grid}
    best = min(profiles,key=profiles.get)
    independent = optimize(gp,0,shared=False)
    return {"width_days":width,"microlensing_length_days":ell,
            "n_bin_observations":len(gp.t),"fixed_hyper_best_lag_days":best,
            "nll_common_at_210":round(float(profiles[210]),4),
            "nll_independent_images":round(float(independent.fun),4),
            "nll_independent_minus_common_at_210":round(float(independent.fun-profiles[210]),4),
            "nll_zero_minus_210":round(float(profiles[0]-profiles[210]),4),
            "nll_365_minus_210":round(float(profiles[365]-profiles[210]),4)}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args=parser.parse_args()
    for width in ((14,) if args.quick else (10,14,20,30)):
        for ell in ((180,800) if args.quick else (180,400,800,1200)):
            print(json.dumps(stress(width,ell),sort_keys=True))
