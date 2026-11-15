"""Joint image-delay Gaussian process, including independent extrinsic processes.

Reproduces the new *sensitivity test*, not a posterior interval. A shares the
same intrinsic OU process as B, evaluated at t - delta. Independent squared-
exponential image terms represent microlensing. Mean magnitudes are integrated
with flat priors (generalized least squares). The independent-images control
deletes the cross-image intrinsic covariance entirely.

Run: python results/j0826-lens-delay/followup/gaussian_process.py [--quick]
Requires numpy, scipy. Uses immutable CSVs already in PR #22.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve, LinAlgError
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
TARGET = "J0826+7002"


def observations(instrument="ztf", width=14):
    if instrument == "ztf":
        with (ROOT / "ztf-nights.csv").open(newline="") as f:
            rows = [r for r in csv.DictReader(f) if r["system"] == TARGET and r["band"] == "zg"]
        out = []
        for image in ("B", "A"):
            blocks = defaultdict(list)
            for r in rows:
                if r["image"] == image:
                    t = float(r["mjd"])
                    blocks[int(t // width)].append((t, float(r["magnitude"]), float(r["median_catalog_error"])))
            for k in sorted(blocks):
                group = np.array(blocks[k])
                t, mag, e = np.median(group, axis=0)
                # Each disjoint time-bin is one datum, not every original exposure.
                e = np.sqrt(e**2 / len(group) + (.02 if image == "A" else .04)**2)
                out.append((t, mag, e, float(image == "A")))
        return np.array(out, float)
    if instrument == "gaia":
        with (ROOT / "gaia-nights.csv").open(newline="") as f:
            rows = [r for r in csv.DictReader(f) if r["band"] == "g"]
        return np.array([(float(r["bjd_tcb_minus_2400000p5"]), float(r["magnitude"]),
                          max(.02, float(r["median_catalog_error"])), float(r["image"] == "A"))
                         for r in rows], float)
    raise ValueError("Unknown instrument")


class DelayedGP:
    def __init__(self, observations):
        self.d = observations
        self.t, self.y, self.e, im = observations.T
        self.im = im.astype(bool)  # A, which arrives later
        self.eq = self.im[:, None] == self.im[None, :]
        self.aa = self.im[:, None] & self.im[None, :]
        self.bb = (~self.im[:, None]) & (~self.im[None, :])
        self.dist = abs(self.t[:, None] - self.t[None, :])
        self.mean_design = np.column_stack([~self.im, self.im]).astype(float)

    def covariance(self, delta, p, shared=True):
        sigma, tau, mu_a, mu_b, ell, floor_a, floor_b = np.exp(p)
        source_time = self.t - delta * self.im
        intrinsic = sigma**2 * np.exp(-abs(source_time[:, None] - source_time[None, :]) / tau)
        if not shared:
            intrinsic *= self.eq
        extrinsic = np.exp(-.5 * (self.dist / ell)**2) * (mu_a**2 * self.aa + mu_b**2 * self.bb)
        K = intrinsic + extrinsic
        K[np.diag_indices_from(K)] += self.e**2 + np.where(self.im, floor_a**2, floor_b**2) + 1e-8
        return K

    def score(self, delta, p, shared=True):
        try:
            K = self.covariance(delta, p, shared=shared)
            ch = cho_factor(K, lower=True, check_finite=False)
            SinvX = cho_solve(ch, self.mean_design, check_finite=False)
            C = self.mean_design.T @ SinvX
            invC = np.linalg.inv(C)
            beta = invC @ (self.mean_design.T @ cho_solve(ch, self.y, check_finite=False))
            residual = self.y - self.mean_design @ beta
            d = residual @ cho_solve(ch, residual, check_finite=False)
            logdet = 2 * np.log(np.diag(ch[0])).sum()
            logdetC = np.linalg.slogdet(C)[1]
            return float(.5 * (d + logdet + logdetC + (len(self.y) - 2) * np.log(2 * np.pi)))
        except (LinAlgError, ValueError, OverflowError):
            return 1e99

    def heldout_conditional_nll(self, split_mjd, delta, p):
        """Proper predictive density for late observations conditional on early ones.
        Fits means from the early period and marginalizes their uncertainty.
        A single common hyperparameter vector is locked from the early period.
        """
        early = np.flatnonzero(self.t < split_mjd)
        late = np.flatnonzero(self.t >= split_mjd)
        K = self.covariance(delta, p)
        Ktr = K[np.ix_(early, early)]
        Xtr = self.mean_design[early]
        ch = cho_factor(Ktr, lower=True, check_finite=False)
        invX = cho_solve(ch, Xtr, check_finite=False)
        invy = cho_solve(ch, self.y[early], check_finite=False)
        Cinv = np.linalg.inv(Xtr.T @ invX)
        beta = Cinv @ (Xtr.T @ invy)
        cross = K[np.ix_(late, early)]
        xt = self.mean_design[late]
        corrected = xt - cross @ invX
        mean = xt @ beta + cross @ cho_solve(ch, self.y[early] - Xtr @ beta, check_finite=False)
        covariance = K[np.ix_(late, late)] - cross @ cho_solve(ch, cross.T, check_finite=False) + corrected @ Cinv @ corrected.T
        check = cho_factor(covariance, lower=True, check_finite=False)
        residual = self.y[late] - mean
        return float(.5 * (residual @ cho_solve(check, residual, check_finite=False)
                    + 2 * np.log(np.diag(check[0])).sum() + len(late) * np.log(2*np.pi)))


INITIAL = np.log([.12, 160., .35, .09, 650., .045, .075])
BOUNDS = np.log(np.array([[.03, .5], [30, 1400], [.04, 1.],
                          [.01, .5], [180, 1600], [.01, .3], [.01, .35]]))


def optimize(gp, delta, shared=True):
    return minimize(lambda p: gp.score(delta, p, shared), INITIAL,
                    method="L-BFGS-B", bounds=BOUNDS, options={"maxiter": 45, "ftol": 2e-6})


def report(instrument, quick=False):
    gp = DelayedGP(observations(instrument))
    candidates = [0, 210, 215, 235, 365] if quick else [-210, 0, 110, 170, 190, 200, 210, 215, 220, 230, 235, 250, 265, 280, 300, 330, 365]
    fitted = {d: optimize(gp, d) for d in candidates}
    ref = fitted[210].x
    grid = np.arange(-100, 601, 5)
    profile = [(int(d), float(gp.score(d, ref))) for d in grid]
    peak = min(profile, key=lambda t: t[1])
    out = {
        "instrument": instrument, "n": len(gp.t), "bins": 14 if instrument == "ztf" else None,
        "best_fixed_hyperparameters": peak, "refit_nll": {str(d): round(r.fun, 5) for d, r in fitted.items()},
        "independent_images_null_nll": round(optimize(gp, 0, shared=False).fun, 5),
        "hyperparameters_fit_at_210": np.exp(ref).round(5).tolist(),
    }
    if instrument == "ztf":
        out["locked_hyperparameter_future_holdouts"] = {}
        for split in (59600, 60000, 60200):
            train = DelayedGP(gp.d[gp.t < split])
            hyp = optimize(train, 210).x
            scores = {str(d): round(gp.heldout_conditional_nll(split, d, hyp), 4)
                      for d in (0, 210, 215, 225, 235, 365, -210)}
            out["locked_hyperparameter_future_holdouts"][str(split)] = scores
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    for instrument in ("ztf", "gaia"):
        print(json.dumps(report(instrument, args.quick), sort_keys=True))


if __name__ == "__main__":
    main()
