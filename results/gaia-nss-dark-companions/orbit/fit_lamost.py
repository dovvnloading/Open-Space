"""Orbit fit conditional on a dark companion for Gaia DR3 3378588057203660160."""
import csv
import hashlib
import json

import numpy as np

from fit_controls import HERE, DR3_JD, fit, nss_constraint

SID = "3378588057203660160"
NSS = HERE/"lamost_candidate_nss.csv"
RV = HERE.parent/"lamost"/"exposures.csv"
MASS = HERE.parent/"gaia"/"mass_screen.csv"


def main():
    with NSS.open() as f:
        row = next(csv.DictReader(f))
        assert row["source_id"]==SID and row["nss_solution_type"]=="Orbital"
    with RV.open() as f:
        exposures = [r for r in csv.DictReader(f) if r["source_id"]==SID]
    with MASS.open() as f:
        mass=next(r for r in csv.DictReader(f) if r["source_id"]==SID)
    m1 = float(mass["m1"])
    sigma = (float(mass["m1_upper"])-float(mass["m1_lower"]))/2
    nights=[]
    for date in sorted({r["obsdate"] for r in exposures}):
        group=[r for r in exposures if r["obsdate"]==date]
        t,y,s=np.array([[float(r[k]) for k in ["mjd_utc_start","rv_kms","err_kms"]] for r in group]).T
        w=s**-2
        nights.append([np.average(t,weights=w)+2400000.5-DR3_JD,
                       np.average(y,weights=w),1/np.sqrt(w.sum())])
    t,y,s=np.array(nights).T
    all_fits={}
    rng=np.random.default_rng(374452)
    with (HERE/"lamost_candidate_residuals.csv").open("w") as f:
        writer=csv.writer(f)
        writer.writerow(["fit","jd_rel_reference","instrument","rv_kms","err_kms","model_kms",
                         "observed_minus_model_kms"])
        for nss_scale, floor, primary in [(1,1,(m1,sigma)),(3,2,(m1,sigma)),(1,1,(1.8141,.0546))]:
            label=f"NSSx{nss_scale}_floor{floor}_M1{primary[0]:.3f}"
            astro,start=nss_constraint(row,nss_scale)
            start[8]=np.median(y)
            all_fits[label]=fit(label,start,astro,t,y,np.hypot(s,floor),["LAMOST"]*len(t),
                                 rng,writer,primary_mass=primary,include_mass_density=True)
    # Fix the temporal split at the start of the final observing season.
    # The last seven nights are never used while fitting this prediction.
    keep=(t+DR3_JD)<2459100.
    astro,start=nss_constraint(row,1.)
    start[8]=np.median(y[keep])
    err=np.hypot(s,1.)
    with (HERE/"lamost_training_residuals.csv").open("w") as f:
        writer=csv.writer(f)
        writer.writerow(["fit","jd_rel_reference","instrument","rv_kms","err_kms","model_kms",
                         "observed_minus_model_kms"])
        split=fit("through_March_2020",start,astro,t[keep],y[keep],err[keep],
                  ["LAMOST"]*sum(keep),rng,writer,primary_mass=(m1,sigma),
                  holdout=(t[~keep],y[~keep],err[~keep]))
    details=dict(holdout_fit=split,source_id=SID,exposures=len(exposures),baseline_d=np.ptp(t),
                 sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [NSS,RV,MASS]},
                 models=all_fits)
    (HERE/"lamost_candidate_fits.json").write_text(json.dumps(details,indent=2)+"\n")
if __name__=="__main__":
    main()
