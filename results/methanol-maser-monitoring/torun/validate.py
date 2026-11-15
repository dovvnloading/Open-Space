"""Evaluate Hitachi-trained periods and maximum epochs on Torun vector photometry."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.stats import t as student_t

HERE = Path(__file__).resolve().parent
FROZEN = {16.13: (112.76950, 59062.28), 16.92: (112.75746, 59067.65)}


def signed_fit(y, X):
    beta = np.linalg.lstsq(X,y,rcond=None)[0]
    dof = len(y)-len(beta)
    err = np.sqrt(np.linalg.inv(X.T @ X)[-1,-1] * np.sum((y-X@beta)**2)/dof)
    return {"amplitude_jy": beta[-1], "se_jy": err,
            "p_one_sided": student_t.sf(beta[-1]/err,dof), "n": len(y)}


def red_noise_test(t, y, baseline, cosine, rng):
    lag = abs(t[:,None]-t[None,:])

    def likelihood(theta, fit=False):
        white, red, tau = np.exp(theta)
        K = np.eye(len(t))*white**2 + red**2*np.exp(-lag/tau)
        cf = cho_factor(K,lower=True)
        KX = cho_solve(cf,baseline)
        coeff = np.linalg.solve(baseline.T @ KX, KX.T @ y)
        res = y - baseline@coeff
        if fit:
            return baseline@coeff, cf, K
        return np.log(np.diag(cf[0])).sum() + .5*res@cho_solve(cf,res)

    fits = [minimize(likelihood, np.log([.3,.4,k]), method="L-BFGS-B",
                     bounds=[(-7,2),(-7,2),(-1,8)]) for k in (1,10,60,180,800)]
    best = min((fit for fit in fits if fit.success), key=lambda fit: fit.fun)
    mu, cf, K = likelihood(best.x, fit=True)
    X = np.column_stack([baseline,cosine])
    KX = cho_solve(cf,X)
    cov = np.linalg.inv(X.T@KX)
    weights = (cov@KX.T)[-1]
    amplitude = weights@y
    null = rng.multivariate_normal(mu,K,20000)@weights
    white, red, tau = np.exp(best.x)
    return {"white_jy":white, "OU_jy":red,"OU_decay_days":tau,"amplitude_jy":amplitude,
            "conditional_se_jy":np.sqrt(cov[-1,-1]),"simulated_p_one_sided":
            (1+int(sum(null>=amplitude)))/(len(null)+1)}


def test(series, period, epoch):
    t = series.mjd.to_numpy()
    y = series.flux_jy.to_numpy()
    n = len(t)
    phase = 2*np.pi*(t-epoch)/period
    baseline = np.column_stack([np.ones(n),(t-t.mean())/np.ptp(t)])
    cosine = np.cos(phase)
    X = np.column_stack([baseline,cosine])
    model = signed_fit(y,X)
    model["halves"] = [signed_fit(y[m],X[m]) for m in (t<np.median(t),t>=np.median(t))]
    rng = np.random.default_rng(20261115)
    model["red_noise"] = red_noise_test(t,y,baseline,cosine,rng)
    residual = lambda a: a - baseline@np.linalg.lstsq(baseline,a,rcond=None)[0]
    r, z = residual(y), residual(cosine)
    model["moving_block_bootstrap"] = {}
    for block in (2,3,4,6):
        starts = rng.integers(n,size=(30000,int(np.ceil(n/block))))
        ids = (starts[:,:,None]+np.arange(block))%n
        null = r[ids.reshape(30000,-1)[:,:n]] @ z/(z@z)
        model["moving_block_bootstrap"][block] = (1+int(sum(null >= model["amplitude_jy"])))/(len(null)+1)
    # Only after the prespecified signed test, estimate a free phase offset.
    free = np.linalg.lstsq(np.column_stack([X,np.sin(phase)]),y,rcond=None)[0]
    model.update({"period_days":period,"Hitachi_maximum_mjd":epoch,
                  "descriptive_phase_shift_days":np.arctan2(free[-1],free[-2])*period/(2*np.pi)})
    return model


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",type=Path,default=HERE/"torun.csv")
    p.add_argument("--output",type=Path,default=HERE/"fixed_phase.json")
    a=p.parse_args()
    groups=pd.read_csv(a.input).sort_values("mjd").groupby("velocity_km_s")
    result={str(v):test(g,*FROZEN[v]) for v,g in groups}
    a.output.write_text(json.dumps(result,indent=2)+"\n")
