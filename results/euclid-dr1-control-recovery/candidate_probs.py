#!/usr/bin/env python3
"""Conditional two-class (known quasar vs labelled dwarf) colour likelihoods.

Not a quasar posterior for a survey field: tuned selection, biased labels,
unmeasured galaxy/blend/artifact components, and a *hypothesized* class prior.
Never extrapolate numeric probabilities beyond empirical labelled support.
"""
import json,math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
BANDWIDTH=0.32
PRIOR=0.001
ALT_PRIOR=0.01
MAX_DISTANCE=4.0
def load(n):return json.loads((OUT/n).read_text())
def save(n,a):(OUT/n).write_text(json.dumps(a,indent=2,sort_keys=True,allow_nan=False)+"\n")
def points():
    c=load("control_audit.json")
    q=[];d=[]
    for r in c:
        if "source" not in r:continue
        y=r["derived"].get("Y_minus_J_AB")
        if y is None and r["source"].get("flux_y_3fwhm_aper",0)<=0:
            y=r["derived"].get("Y_minus_J_2sig_limit")
        h=r["derived"].get("J_minus_H_AB")
        if y is None or h is None or not all(math.isfinite(float(x)) for x in (y,h)):
            continue
        (q if r["class"]=="quasar" else d).append([float(y),float(h)])
    if (len(q),len(d))!=(24,14):raise AssertionError("Changed matched labelled control sample")
    return np.array(q),np.array(d)
def density(x,points):
    return float(np.mean(np.exp(-0.5*np.sum(((points-np.asarray(x))/BANDWIDTH)**2,axis=1))))/(2*np.pi*BANDWIDTH**2)
def posterior(x,q,d,p):
    lq,ld=density(x,q),density(x,d)
    den=p*lq+(1-p)*ld
    return p*lq/den if den>0 else None
def annotate(row,q,d):
    y=row.get("Y_minus_J_used")
    h=row.get("J_minus_H_AB")
    out=dict(row)
    out.update(conditional_p_quasar_vs_dwarf=None,conditional_p_dwarf_vs_quasar=None,
               conditional_p_quasar_prior_1pct=None,
               conditional_colour_density_status="UNAVAILABLE",
               color_distance_from_label_in_bandwidths=None,
               field_quasar_probability=None)
    if y is None or h is None:return out
    x=np.array([float(y),float(h)])
    if not np.isfinite(x).all():return out
    dist=float(np.sqrt(np.sum(((np.concatenate((q,d))-x)/BANDWIDTH)**2,axis=1)).min())
    out["color_distance_from_label_in_bandwidths"]=round(dist,4)
    if dist>MAX_DISTANCE:
        out["conditional_colour_density_status"]="OUTSIDE_EMPIRICAL_LABELLED_SUPPORT_NO_PROBABILITY"
        return out
    a=posterior(x,q,d,PRIOR);b=posterior(x,q,d,ALT_PRIOR)
    if a is None or b is None:
        out["conditional_colour_density_status"]="ZERO_DENSITY_UNDEFINED"
        return out
    out["conditional_colour_density_status"]="UNVALIDATED_CONDITIONAL_TWO_CLASS_MODEL"
    out["conditional_p_quasar_vs_dwarf"]=round(a,8)
    out["conditional_p_dwarf_vs_quasar"]=round(1-a,8)
    out["conditional_p_quasar_prior_1pct"]=round(b,8)
    return out
def main():
    q,d=points()
    existing=load("top_candidates_with_checks.json")
    scored=[annotate(r,q,d) for r in existing]
    save("ranked_conditional_top100.json",scored)
    supported=sorted((r for r in scored if r["conditional_p_quasar_vs_dwarf"] is not None),
                     key=lambda r:(-r["conditional_p_quasar_vs_dwarf"],
                                   -r["rank_score"],r["id"]))
    save("conditional_supported_top100.json",supported)
    loo=[]
    for group,points,other,truth in (("quasar",q,d,1),("dwarf",d,q,0)):
        for i,x in enumerate(points):
            same=np.delete(points,i,axis=0)
            proba=posterior(x,same,other,PRIOR) if truth else posterior(x,other,same,PRIOR)
            if proba is None:raise AssertionError("Missing LOO probability")
            loo.append({"class":group,"p":proba,"truth":truth})
    brier=sum((t["p"]-t["truth"])**2 for t in loo)/len(loo)
    model={
        "status":"UNVALIDATED_Q_VS_DWARF_ONLY",
        "prior_Q_given_Q_or_dwarf":PRIOR,
        "alternate_prior_Q_given_Q_or_dwarf":ALT_PRIOR,
        "kernel":"bivariate isotropic Gaussian KDE, equal weights within each class",
        "bandwidth_AB_mag":BANDWIDTH,
        "numeric_rejected_when_nearest_control_gt_bandwidths":MAX_DISTANCE,
        "training_label_quasars":len(q),"training_label_dwarfs":len(d),
        "cross_validated_leave_one_out_Brier":round(brier,6),
        "top100_with_numeric_conditional_probability":len(supported),
        "top100_outside_colour_support":sum(
           r["conditional_colour_density_status"]=="OUTSIDE_EMPIRICAL_LABELLED_SUPPORT_NO_PROBABILITY"
           for r in scored),
        "conditional_labels_only":True,
        "field_quasar_probability_calibrated":False,
        "field_completeness_calibrated":False,
        "selection_tuned_using_same_controls":True,
        "not_modelled":["galaxies","blends","artifacts","other stars","survey selection function","source coverage","spectral confirmation"],
        "calibration_limitation":"No field purity, posterior, or usable discovery-significance inference follows from this two-class labelled KDE.",
        "new_spectroscopic_confirmations":0,
    }
    save("conditional_model.json",model)
    print("CONDITIONAL_MODEL",json.dumps(model),flush=True)
if __name__=="__main__":main()
