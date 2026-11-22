"""Post-selection, paired-reduction diagnostics for the 72 d HD 17970 signal."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from runtime import require_survey_workflow
from search import nightly


STAR = "GJ3187"
YEAR = 365.25
JULY_2000 = 2451726.5
T0 = 2455000.0  # Fixed time reference for reported amplitudes and phases
TAU = 2 * np.pi


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def exposures(path, name):
    chunks = []
    for df in pd.read_csv(path, chunksize=25000, dtype={"star": str, "epoch_id": str}):
        chunks.append(df[df.star.eq(name)])
    return pd.concat(chunks, ignore_index=True).sort_values("bjd").reset_index(drop=True)


def matched(a, b):
    delta = 2 / 86400
    x, y = a.bjd.to_numpy(), b.bjd.to_numpy()
    def slots(v, w):
        return np.searchsorted(w, v-delta), np.searchsorted(w, v+delta, side="right")
    l, h = slots(x, y)
    lb, hb = slots(y, x)
    i = np.flatnonzero(h-l == 1)
    j = l[i]
    keep = (hb[j]-lb[j] == 1) & (lb[j] == i)
    aa, bb = a.iloc[i[keep]].reset_index(drop=True), b.iloc[j[keep]].reset_index(drop=True)
    valid = aa.era.eq(bb.era) & np.floor(aa.bjd).eq(np.floor(bb.bjd))
    identifiers = aa.epoch_id.notna() & bb.epoch_id.notna() & ~aa.epoch_id.duplicated(False) & ~bb.epoch_id.duplicated(False)
    return aa[valid & identifiers], bb[valid & identifiers]


def design(t, era, season, train, periods, base):
    dt = t - T0
    cols = [dt / YEAR]
    levels = era if base == "era_trend" else np.array([f"{s}/{e}" for s,e in zip(season, era)])
    for level in np.unique(levels[train]):
        cols.append((levels == level).astype(float))
    for period in periods:
        cols.extend([np.cos(TAU*dt/period), np.sin(TAU*dt/period)])
    return np.column_stack(cols)


def solve(x, y, err, train):
    a = x[train] / err[train,None]
    u, s, vt = np.linalg.svd(a, full_matrices=False)
    if not len(s) or np.sum(train) <= x.shape[1] + 1 or s[-1] <= s[0]*1e-10:
        return None
    pinv = (vt.T / s) @ u.T
    beta = pinv @ (y[train] / err[train])
    resid = (y - x @ beta) / err
    rss = np.sum(resid[train]**2)
    cov = (pinv @ pinv.T) * max(1., rss / (sum(train)-len(beta)))
    return beta, resid, rss, cov, s[0]/s[-1]


def inspect_fluxes(obs, source, models, periods, model_rows, folds, components):
    t, y, err, era = (obs[c].to_numpy() for c in ["bjd", "rv", "err", "era"])
    season = np.floor((t-JULY_2000)/YEAR).astype(int)
    all_obs = np.ones(len(t), dtype=bool)
    for base in ["era_trend","season_era_trend"]:
        for label, (names, parent) in models.items():
            selected = [periods[n] for n in names]
            row = dict(reduction=source, nuisance=base, model=label, compare_to=parent,
                       n=len(t), score_kind="absolute" if base=="era_trend" else "within-season-era contrasts")
            model_rows.append(row)
            x = design(t, era, season, all_obs, selected, base)
            fit = solve(x,y,err,all_obs)
            if fit is None:
                row["outcome"] = "insufficient rank or nights"
                continue
            beta, residual, rss, cov, cond = fit
            row.update(outcome="fit", chi2=rss, n_parameters=len(beta), bic=rss+len(beta)*np.log(len(t)),
                       design_condition=cond)
            for k, name in enumerate(names):
                start = len(beta)-2*len(names)+2*k
                a,b = beta[start:start+2]
                v = cov[start:start+2,start:start+2]
                amp = np.hypot(a,b)
                cc = dict(reduction=source,nuisance=base,model=label,term=name,period_d=periods[name],
                          ref_bjd=T0,cos_ms=a,sin_ms=b,cos_var=v[0,0],sin_var=v[1,1],cos_sin_cov=v[0,1],k_ms=amp)
                if amp>0:
                    grad_amp = np.array([a,b])/amp
                    grad_phase = np.array([-b,a])/(amp**2*TAU)
                    cc.update(k_se_ms=np.sqrt(max(0,grad_amp@v@grad_amp)),
                              max_phase_cycles=np.arctan2(b,a)/TAU % 1,
                              phase_se_cycles=np.sqrt(max(0,grad_phase@v@grad_phase)))
                components.append(cc)
            for s in np.unique(season):
                test = season == s
                train = ~test
                frow = dict(reduction=source,nuisance=base,model=label,season=int(s),n_test=int(test.sum()))
                folds.append(frow)
                if base == "era_trend" and not set(era[test]).issubset(era[train]):
                    frow["outcome"] = "untrained fibre era"
                    continue
                xx = design(t,era,season,train,selected,base)
                trained = solve(xx,y,err,train)
                if trained is None:
                    frow["outcome"] = "insufficient training rank or nights"
                    continue
                # No held-out velocity enters training. For the seasonal model
                # the held-out observables are contrasts orthogonal to an unknown
                # season/era intercept, eliminating one degree per such group.
                prediction = xx[test] @ trained[0]
                r = y[test]-prediction
                dof = len(r)
                if base == "season_era_trend":
                    for e in np.unique(era[test]):
                        g = era[test]==e
                        r[g] -= np.average(r[g], weights=err[test][g]**-2)
                        dof -= 1
                if dof == 0:
                    frow["outcome"] = "no within-season contrast"
                    continue
                frow.update(outcome="scored",score_dof=dof,
                            predictive_chi2=np.sum((r/err[test])**2))
    # Period and sample were chosen with all discovery data. Scores diagnose
    # stability and never become a fresh FAP or independent confirmation.


def compare(rows, folds):
    frame = pd.DataFrame(rows)
    cv = pd.DataFrame(folds)
    for i,row in frame.iterrows():
        ref = frame[(frame.reduction==row.reduction)&(frame.nuisance==row.nuisance)&(frame.model==row.compare_to)]
        if ref.empty or ref.iloc[0].outcome!="fit" or row.outcome!="fit":
            continue
        ref = ref.iloc[0]
        frame.loc[i,"delta_bic_vs_parent"] = ref.bic - row.bic
        frame.loc[i,"conditional_chi2_gain"] = ref.chi2 - row.chi2
        frame.loc[i,"conditional_fraction"] = (ref.chi2-row.chi2)/ref.chi2 if ref.chi2>0 else np.nan
        same = cv[(cv.reduction==row.reduction)&(cv.nuisance==row.nuisance)&(cv.outcome=="scored")]
        pair = same[same.model==row.model].merge(same[same.model==row.compare_to],on="season")
        frame.loc[i,"common_scored_seasons"] = len(pair)
        if len(pair):
            gain = pair.predictive_chi2_y-pair.predictive_chi2_x
            frame.loc[i,"common_score_dof"] = pair.score_dof_x.sum()
            frame.loc[i,"cv_chi2_gain_vs_parent"] = gain.sum()
            frame.loc[i,"cv_seasons_improved"] = sum(gain>0)
    return frame


def run(args, receipt):
    h = pd.read_csv(args.evidence / "holdout.csv",keep_default_na=False)
    v = pd.read_csv(args.evidence / "signals.csv")
    row = h[h.star.eq(STAR) & h.signal_index.eq(1)]
    if len(row)!=1:
        raise ValueError("unique target missing from holdout roster")
    target = row.iloc[0]
    if (str(target.verified_simbad_oid).split(".")[0]!="1188552" or
            target.identity_position!="vetted SIMBAD J2000" or
            "".join(target.eso_id.split())!="HD17970" or not 0<=float(target.sep_arcsec)<=3):
        raise ValueError("archived reciprocal name and coordinate check failed")
    receipt.update(eso_identity=target.eso_id,simbad_oid=1188552,identity_sep_arcsec=float(target.sep_arcsec))
    fits = v[v.star.eq(STAR)].set_index("signal_index",verify_integrity=True)
    pc,pf = fits.loc[[1,2],"period_d"]
    if not (72<pc<73 and 360<pf<365):
        raise ValueError("discovery and comparison periods differ from proposal")
    p = {"candidate":pc,"other":pf,"other5":pf/5,**{f"year{k}":YEAR/k for k in range(1,6)}}
    periods = list(p)[3:]
    models = {"null":([],""),"candidate":(["candidate"],"null"),"other":(["other"],"null"),
              "other5":(["other5"],"null"),"other_and_fifth":(["other","other5"],"other")}
    for term in periods:
        models[term]=([term],"null")
    for name, nuisance in [("year5",["year5"]),("other5",["other5"]),
                           ("other_and_fifth",["other","other5"]),("annual_all",periods),
                           ("annual_all_other_fifth",periods+["other","other5"])]:
        if name not in models:
            models[name]=(nuisance,"annual_all" if name=="annual_all_other_fifth" else "null")
        models[name+"_candidate"]=(nuisance+["candidate"],name)
    a, b = exposures(args.primary,STAR), exposures(args.eso,target.eso_id)
    aa,bb = matched(a,b)
    receipt.update(primary_epochs=len(a),eso_epochs=len(b),matched_exposures=len(aa),periods_d=p)
    ids = pd.concat([aa[["epoch_id","bjd","rv","err","era"]].reset_index(drop=True).add_prefix("bank_"),
                     bb[["epoch_id","bjd","rv","err","era"]].reset_index(drop=True).add_prefix("eso_")],axis=1)
    ids.to_csv(args.out/"deep-paired-exposures.csv",index=False)
    if len(aa)<12:
        raise ValueError("fewer than twelve uniquely matched exposures")
    na, nb = nightly(aa),nightly(bb)
    if not na[["day","era"]].equals(nb[["day","era"]]) or len(na)<12:
        raise ValueError("night pairing failed or fewer than twelve paired nights")
    receipt.update(paired_nights=len(na))
    paired = na[["day","era","bjd","rv","err"]].merge(nb[["day","era","bjd","rv","err"]],on=["day","era"],
                                                      suffixes=("_bank","_eso"),validate="one_to_one")
    paired.to_csv(args.out/"deep-paired-nights.csv",index=False)
    rows,folds,components=[],[],[]
    for obs,src in [(na,"RVBank SERVAL"),(nb,"ESO DRS")]:
        inspect_fluxes(obs,src,models,p,rows,folds,components)
    compare(rows,folds).to_csv(args.out/"deep-models.csv",index=False)
    pd.DataFrame(folds).to_csv(args.out/"deep-season-scores.csv",index=False)
    pd.DataFrame(components).to_csv(args.out/"deep-components.csv",index=False)


def main():
    require_survey_workflow()
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ["primary","eso","evidence","out"]:
        parser.add_argument("--"+option,required=True,type=Path)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    files=[args.primary,args.eso,args.evidence/"holdout.csv",args.evidence/"signals.csv",Path(__file__)]
    receipt=dict(star=STAR,inputs={str(f):sha(f) for f in files if f.exists()},method=(
        "Reciprocal unique exposure within two seconds, same fibre era and integer BJD night. "
        "Standard independent nightly weights with internal-scatter check and 1 m/s floor. "
        "Fixed periods, weighted linear sine/cosine and common slope. Months grouped into July–June observing years. "
        "An intercept per era or per observing year/era. BIC assumes the stated independent Gaussian errors. "
        "Positive conditional chi-square, BIC and predictive gains favour the larger model. Leave-one-season-out "
        "uses only training velocities for coefficients. Seasonal scores concern within-season/era contrasts. "
        "Phase is time of cosine maximum after BJD 2455000 in cycles. Formal first-order covariance conditional "
        "on fixed periods, scaled by max(1,reduced chi2). Weak signals and nearly collinear harmonics have unreliable "
        "phase linearizations. Same spectra, post-selection epochs and periods, possible shared systematics and "
        "unmodelled correlated noise preclude an independent detection or fresh false-alarm probability."))
    try:
        run(args,receipt)
        receipt["status"]="completed"
    except (ValueError,OSError,KeyError,np.linalg.LinAlgError) as exc:
        receipt.update(status="unavailable",reason=str(exc))
    (args.out/"deep-followup.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(receipt["status"],receipt.get("reason",""))


if __name__=="__main__":
    main()
