"""Compare frozen HARPS predictions with correlated ESPRESSO velocity residuals."""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent


def fit_noise(t, velocity, errors, model):
    """Profile an offset and linear drift, then fit the noise hyperparameters."""
    dt = np.abs(t[:, None] - t)
    design = np.column_stack((np.ones(len(t)), (t - t.mean()) / 365.25))
    bounds = [(np.log(.05), np.log(15.))]
    starts = [[2.0]]
    if model != "white":
        bounds += [(np.log(.05), np.log(15.)), (np.log(2.), np.log(300.))]
        starts = [[2., 3., lag] for lag in (5., 20., 90.)]
    if model == "quasiperiodic":
        bounds += [(np.log(15.), np.log(65.))]
        starts = [row + [rotation] for row in starts for rotation in (20., 24., 32., 40., 49., 57.)]

    def objective(pars):
        values = np.exp(pars)
        cov = np.diag(errors**2 + values[0]**2)
        if model == "exponential":
            cov += values[1]**2 * np.exp(-dt/values[2])
        if model == "quasiperiodic":
            # Fixed phase smoothness of one radian. The stellar period and
            # decorrelation time are allowed to vary under each hypothesis.
            cov += values[1]**2 * np.exp(-.5*(dt/values[2])**2
                     - 2*np.sin(np.pi*dt/values[3])**2)
        factor = cho_factor(cov, lower=True, check_finite=False)
        ci_x = cho_solve(factor, design, check_finite=False)
        coeff = np.linalg.solve(design.T @ ci_x, ci_x.T @ velocity)
        whitened = solve_triangular(factor[0], velocity-design @ coeff,
                                    lower=True, check_finite=False)
        return float(np.log(np.diag(factor[0])).sum() + .5 * (whitened @ whitened))

    fits = [minimize(objective, np.log(init), method="L-BFGS-B", bounds=bounds,
                     options={"maxiter": 150}) for init in starts]
    fits = [f for f in fits if f.success]
    best = min(fits, key=lambda f: f.fun)
    return best.fun, np.round(np.exp(best.x), 3).tolist()


def main():
    with (ROOT / "espresso.csv").open() as f:
        rows = list(csv.DictReader(f))
    t = np.asarray([float(r["bjd"])-2456000 for r in rows])
    y = np.asarray([1000*float(r["rv"]) for r in rows])
    y -= y.mean()
    e = np.asarray([1000*float(r["rv_err"]) for r in rows])
    models = json.loads((ROOT / "harps_predictions.json").read_text())
    best_harps = min(models, key=lambda m: m["rms"])
    phase = 2*np.pi*t[:, None]/best_harps["P"]
    coefficients = np.asarray(best_harps["params"][:8]).reshape(4,2)
    terms = np.cos(phase)*coefficients[:,0] + np.sin(phase)*coefficients[:,1]
    choices = {"null": [], "2.56_only": [0], "both_inner": [0,1],
               "four": [0,1,2,3], "omit_2.56": [1,2,3], "omit_11.06": [0,2,3]}
    output = {}
    for model in ("white", "exponential", "quasiperiodic"):
        result = {}
        for name, use in choices.items():
            score, pars = fit_noise(t, y-terms[:,use].sum(axis=1), e, model)
            result[name] = (score, pars)
        output[model] = {name: {"delta_log_likelihood_vs_null": round(result["null"][0]-score, 2),
                                 "noise_parameters": pars} for name, (score, pars) in result.items()}
    return output


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
