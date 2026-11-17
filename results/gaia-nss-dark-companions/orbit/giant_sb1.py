"""RV confirmation and primary-dependent edge-on mass of DR3 275582835635329024."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

from fit_controls import HERE, DR3_JD, companion
from sdss_sb1 import FIELDS, GAIA_MIN, GAIA_MAX, covariance, curve

SID="275582835635329024"
EXPOSURES=HERE.parent/"extension"/"lamost"/"exposures.csv"
INPUT=HERE/"flame_giant_nss.csv"
MASS=HERE.parent/"extension"/"mass_screen.csv"


def weighted_q(a,w):
    order=np.argsort(a)
    wc=np.cumsum(w[order])-.5*w[order]
    return np.interp([.16,.50,.84],wc,a[order]).tolist()


def analyze(prior, mu, L, t, y, s, rng, train):
    # Integrate over a uniform additive LAMOST-Gaia zeropoint analytically.
    rv=curve(prior,t)
    offsets=np.average(y[train]-rv[:,train],axis=1,weights=s[train]**-2)
    rv+=offsets[:,None]
    logw=-.5*np.sum(((y[train]-rv[:,train])/s[train])**2,axis=1)
    w=np.exp(logw-logw.max())
    w/=w.sum()
    idx=rng.choice(len(w),size=30000,p=w)
    err_zp=1/np.sqrt(np.sum(s[train]**-2))
    predictions=rv[idx,:]+rng.normal(0,err_zp,size=(len(idx),1))
    def residual(x):
        g=solve_triangular(L,x[:6]-mu,lower=True)
        ext=(curve(x[:6],t).ravel()+x[6]-y)/s
        return np.r_[g,ext[train]]
    fit=least_squares(residual,np.r_[mu,np.median(offsets)],x_scale="jac",
                      bounds=(np.r_[GAIA_MIN,-np.inf],np.r_[GAIA_MAX,np.inf]))
    assert fit.success
    x=fit.x
    fitrv=curve(x[:6],t).ravel()+x[6]
    dof=int(train.sum()-1)
    info=dict(n_training=int(train.sum()),importance_ess=1/(w@w),
              parameter_q16_q50_q84={name:weighted_q(prior[:,i],w) for i,name in enumerate(FIELDS)},
              lamost_offset_kms=weighted_q(offsets+rng.normal(size=len(w))*err_zp,w),
              chi2_gaia=float(np.sum(residual(x)[:6]**2)),chi2_rv=float(np.sum(residual(x)[6:]**2)),
              joint_dof=dof,rv_rms_kms=float(np.sqrt(np.mean((fitrv[train]-y[train])**2))),
              best=dict(zip(FIELDS+["lamost_offset_kms"],x)))
    if not np.all(train):
        test=~train
        sim=predictions[:,test]+rng.normal(size=(len(idx),test.sum()))*s[test]
        mean=sim.mean(axis=0)
        Ltest=np.linalg.cholesky(np.cov(sim.T))
        d=solve_triangular(Ltest,y[test]-mean,lower=True)
        dsim=solve_triangular(Ltest,(sim-mean).T,lower=True)
        info["holdout"]=dict(n=int(test.sum()),predicted_mean_kms=mean.tolist(),
            observed_kms=y[test].tolist(),rms_kms=np.sqrt(np.mean((y[test]-mean)**2)),
            predictive_tail_fraction=(1+(np.sum(dsim**2,axis=0)>=d@d).sum())/(1+len(idx)))
    return info,idx,fitrv


def main(exposures):
    with INPUT.open() as f:r=next(csv.DictReader(f))
    assert r["source_id"]==SID
    with MASS.open() as f:m=next(r for r in csv.DictReader(f) if r["source_id"]==SID)
    with exposures.open() as f:rows=[r for r in csv.DictReader(f) if r["source_id"]==SID]
    dates=sorted({r["obsdate"] for r in rows})
    nights=[]
    for date in dates:
        group=[r for r in rows if r["obsdate"]==date]
        t,y,s=np.array([[float(r[k]) for k in ["mjd_utc_start","rv_kms","err_kms"]] for r in group]).T
        w=s**-2
        nights.append([np.average(t,weights=w)+2400000.5-DR3_JD,np.average(y,weights=w),1/np.sqrt(w.sum())])
    t,y,raw_s=np.array(nights).T
    rng=np.random.default_rng(402768)
    mu=np.array([float(r[k]) for k in FIELDS])
    L=np.linalg.cholesky(covariance(r))
    std=rng.normal(size=(65536,6))
    results={}
    with (HERE/"giant_sb1_residuals.csv").open("w") as f:
        writer=csv.writer(f)
        writer.writerow(["model","date","jd_rel_DR3","rv_kms","err_kms","fitted_rv_kms","residual_kms","in_training"])
        for label,inflate,floor,train in [("nominal",1,1,np.ones(len(t),bool)),
               ("inflated",3,2,np.ones(len(t),bool)),
               ("holdout_last_season",1,1,np.array([d<"2019-10-01" for d in dates]))]:
            s=np.hypot(raw_s,floor)
            prior=mu+std@(L.T*inflate)
            prior=prior[np.all((prior>GAIA_MIN)&(prior<GAIA_MAX),axis=1)]
            info,idx,rv=analyze(prior,mu,L*inflate,t,y,s,rng,train)
            p=prior[idx]
            spec_f=1.036149e-7*p[:,0]*p[:,2]**3*(1-p[:,3]**2)**1.5
            info["mass_function_msun"]=np.quantile(spec_f,[.16,.5,.84]).tolist()
            info["critical_primary_msun_for_min_m2_3"]=np.quantile(np.sqrt(27/spec_f)-3,[.16,.5,.84]).tolist()
            info["mass_models"]={}
            m1=float(m["m1"])
            sig=(float(m["m1_upper"])-float(m["m1_lower"]))/2
            for tag,mean,sigma in [("catalog",m1,sig),("30pct_uncertainty",m1,.30*m1),
                                   ("30pct_lower_primary",.7*m1,0.)]:
                pri=rng.normal(mean,sigma,len(spec_f))
                good=pri>0.
                masses=companion(spec_f[good],pri[good])
                pdf,edges=np.histogram(masses,bins=60,density=True)
                info["mass_models"][tag]=dict(primary_mean=mean,primary_sigma=sigma,
                    density=dict(edges=edges.tolist(),pdf=pdf.tolist()),
                    m2_min_q16_q50_q84=np.quantile(masses,[.16,.5,.84]).tolist(),
                    p_m2_min_gt_1p4=float(np.mean(masses>1.4)),p_m2_min_gt_3=float(np.mean(masses>3)))
            for values in zip(dates,t,y,s,rv,y-rv,train.astype(int)):
                writer.writerow([label,*values])
            results[label]=info
    out=dict(source_id=SID,n_exposures=len(rows),n_nights=len(dates),baseline_d=np.ptp(t),
             input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [INPUT,exposures,MASS]},
             results=results)
    (HERE/"giant_sb1_fits.json").write_text(json.dumps(out,indent=2)+"\n")
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--exposures",type=Path,default=EXPOSURES)
    main(parser.parse_args().exposures)
