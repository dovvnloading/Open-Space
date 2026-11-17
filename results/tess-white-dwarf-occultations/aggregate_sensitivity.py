"""Secondary pointwise upper limit using an aggregate MC bound and zero compatible hosts.

The primary inference in recovery.py is kept separately. This calculation uses
only the 140-source extension selection; the earlier Sector 109 pilot is excluded.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np

def estimate(rows,geometric):
    weights=np.array([float(r["geometry_probability"]) if geometric else 1. for r in rows])
    n=np.array([int(r["trials"]) for r in rows]);k=np.array([int(r["recovered"]) for r in rows])
    expected=float(np.sum(weights*k/n))
    # Each independent uniform [phase, impact] trial adds a bounded term of range w_i/n_i.
    penalty=float(np.sqrt(.5*np.log(1/.005)*np.sum(weights**2/n)))
    lower=max(0.,expected-penalty)
    # AM-GM bounds the probability of zero positives by (1-f*sum_i(q_i)/N)^N.
    upper=len(n)*(1-.045**(1/len(n)))/lower if lower>0 else np.inf
    return dict(effective_stars=expected,mc_lower=lower,mc_penalty=penalty,upper95=float(upper) if upper<1 else None)
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("stars",type=Path);p.add_argument("output",type=Path);a=p.parse_args()
    rows=list(csv.DictReader(a.stars.open()))
    result=[]
    for period,radius in sorted({(float(r["period_days"]),float(r["radius_earth"])) for r in rows}):
        group=[r for r in rows if float(r["period_days"])==period and float(r["radius_earth"])==radius]
        if len(set(r["target"] for r in group))!=len(group):raise ValueError("duplicate targets")
        possible=sum(int(r["possible_host"]) for r in group)
        res=dict(period=period,radius_earth=radius,targets=len(group),possible=possible)
        if possible==0:res.update(transiting=estimate(group,False),intrinsic=estimate(group,True))
        result.append(res)
    a.output.write_text(json.dumps(result,indent=1)+"\n")
