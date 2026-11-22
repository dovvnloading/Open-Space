"""Compare fast weighted periodograms and red covariance to explicit calculations."""
import json
from pathlib import Path
import sys

import numpy as np

from search import Continuum, Search, bounded_uncertainties, fit, survey_trials


def main(output):
    rng=np.random.default_rng(9339141)
    t=np.sort(rng.choice(600, size=420, replace=False)).astype(float) + rng.normal(0,.07,420)
    e=rng.uniform(.6,1.3,t.size)
    y=rng.normal(0,e) + 2.1*np.sin(2*np.pi*t/17.31)
    engine=Search("synthetic-qc","bat",(t,y,e),{},1056)
    powers=engine.raw(y)[0]
    chosen=np.unique(np.r_[np.linspace(0,len(powers)-1,30).astype(int), np.argmax(powers)])
    slow=np.array([fit(engine.t,y,engine.w,engine.freq[k])[1] for k in chosen])
    err=float(np.max(abs(slow-powers[chosen]))/slow.max())
    noise=Continuum(t,y,e,engine.w)
    noise.variances[:]=np.linspace(.01,.12,len(noise.variances))
    noise.configure()
    distance=abs(t[:,None]-t)
    dense=np.diag(noise.white)
    for tau,var in zip(noise.timescales[1:],noise.variances[1:]):
        dense += var*np.exp(-distance/tau)
    project=rng.normal(size=(len(t),4))
    exact=project.T@dense@project
    cov_err=float(np.linalg.norm(noise.bilinear(project)-exact)/np.linalg.norm(exact))
    extreme=e.copy()
    extreme[0]=1e-12
    _,leverage=bounded_uncertainties(t,extreme)
    out=dict(formal_effective_days_with_one_extreme_error=leverage["original_effective_samples"],
             bounded_effective_days_with_one_extreme_error=leverage["effective_samples"],
             two_pass_trials_at_100_frequencies=survey_trials(100,1056),fast_gls_max_difference_over_peak=err,covariance_projection_relative_error=cov_err,
             recovered_period_days=float(1/engine.freq[powers.argmax()]), planted_period_days=17.31)
    print(json.dumps(out,indent=2))
    output.write_text(json.dumps(out,indent=2)+"\n")
    if not (leverage["original_effective_samples"] < 2 and leverage["effective_samples"] > 20 and survey_trials(100,1056)==211200):
        raise RuntimeError("The high-leverage or two-pass trial control failed")
    if err>5e-4 or cov_err>1e-10 or abs(out["recovered_period_days"]-17.31)>.4:
        raise RuntimeError("Fast/exact identity or injected sinusoid calibration failed")


if __name__=="__main__":
    main(Path(sys.argv[1]))
