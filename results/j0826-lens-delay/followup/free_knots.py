"""Joint cubic-spline time-delay model with slow microlensing and movable knots.

Unlike original correlation search, no independent high-pass filtering is done.
One intrinsic spline is evaluated at t_B and t_A-delta; the A-only extrinsic
spline has separate slow knots. B has its own offset. All delays use equal
degrees of freedom and fixed roughness penalties. In the free-knot test,
intrinsic nodes vary by up to 22% of nominal spacing (max 50 days) at *every*
candidate delay. These are sensitivity scores, not calibrated likelihoods.
"""
import argparse
import json
import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import solve
from scipy.optimize import minimize
from gaussian_process import observations

DATA = observations("ztf", 14)
T, Y, E, IM = DATA.T
IM = IM.astype(bool)
SOURCE_RANGE = (58000., 61200.)
MICRO_RANGE = (58000., 61200.)


def knots(step, bounds):
    a, b = bounds
    return np.r_[np.repeat(a, 4), np.arange(a + step, b, step), np.repeat(b, 4)]


def basis(t, k):
    return BSpline.design_matrix(np.clip(t, k[0], k[-1]), k, 3, extrapolate=False).toarray()


def objective(delay, intrinsic, micro, lambda_q=10., lambda_micro=80.):
    q = basis(T - delay*IM, intrinsic)
    m = basis(T, micro)
    nq, nm = q.shape[1], m.shape[1]
    X = np.zeros((len(T), nq + nm + 2))
    X[:, :nq] = q
    X[:, nq:nq+nm] = m*IM[:, None]
    X[:, nq+nm] = ~IM
    X[:, nq+nm+1] = IM
    y = Y - np.where(IM, 17.5, 19.5)
    w = 1/np.sqrt(E**2 + np.where(IM, .045**2, .075**2))
    Z = X*w[:, None]
    R = np.eye(X.shape[1])*1e-6
    R[:nq, :nq] += np.eye(nq)
    R[nq:nq+nm, nq:nq+nm] += np.eye(nm)
    for size, offset, penalty in ((nq, 0, lambda_q), (nm, nq, lambda_micro)):
        for i in range(1, size-1):
            v = np.zeros(X.shape[1])
            v[offset+i-1:offset+i+2] = (1., -2., 1.)
            R += penalty * np.outer(v, v)
    beta = solve(Z.T@Z + R, Z.T@(y*w), assume_a="pos", check_finite=False)
    residual = y*w - Z@beta
    return float(residual@residual + beta@R@beta)


def fixed_grid(qstep=140, mstep=700):
    q, m = knots(qstep, SOURCE_RANGE), knots(mstep, MICRO_RANGE)
    grid = list(range(-80, 561, 5))
    losses = {int(d): objective(d, q, m) for d in grid}
    best = min(losses, key=losses.get)
    return {"intrinsic_knot_days": qstep, "extrinsic_knot_days": mstep, "best_days": best,
            "score_at_zero": round(losses[0], 4), "score_at_210": round(losses[210], 4),
            "score_at_365": round(losses[365], 4), "score_at_best": round(losses[best], 4)}


def free_knots(delay, qstep=190, mstep=700):
    q, m = knots(qstep, SOURCE_RANGE), knots(mstep, MICRO_RANGE)
    ref = q[4:-4].copy()
    allowance = min(50, .22*qstep)
    if not len(ref):
        return objective(delay, q, m)
    def loss(positions):
        moving = q.copy()
        moving[4:-4] = positions
        return objective(delay, moving, m) + .4*sum(((positions-ref)/allowance)**2)
    fit = minimize(loss, ref, method="L-BFGS-B",
                   bounds=[(v-allowance, v+allowance) for v in ref],
                   options={"maxiter":35, "ftol":1e-6})
    return float(fit.fun)


def main():
    arg = argparse.ArgumentParser()
    arg.add_argument("--quick", action="store_true")
    args = arg.parse_args()
    settings = [(90,450),(140,700)] if args.quick else [
        (q,m) for q in (90,140,200,280) for m in (450,700,1000)]
    for q, m in settings:
        print("FIXED", json.dumps(fixed_grid(q,m)))
    trials = [(210,190),(235,190)] if args.quick else [
        (d,q) for q in (130,160,190,240) for d in (0,210,220,235,250,265,280,300,365)]
    for delta, q in trials:
        print("FREE",json.dumps({"delay_days":delta,"intrinsic_spacing_days":q,
                                 "score":round(free_knots(delta,q),5)}))


if __name__ == "__main__":
    main()
