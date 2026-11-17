"""Conditional fit of Gaia's SB1 model and four DR1 visits of a long-period variable."""
import csv
import hashlib
import json
import numpy as np

from fit_controls import HERE,DR3_JD,companion
from sdss_sb1 import FIELDS,GAIA_MIN,GAIA_MAX,covariance
from giant_sb1 import analyze

SID="2951241529971223808"

def main():
    archive=HERE/"desi_lpv_nss.csv"
    visits=HERE.parent/"extension"/"epochs"/"desi_visits.csv"
    masses=HERE.parent/"extension"/"mass_screen.csv"
    with archive.open() as f:r=next(csv.DictReader(f))
    with visits.open() as f:rows=[r for r in csv.DictReader(f) if r["source_id"]==SID and int(r["quality"])]
    with masses.open() as f:m=next(r for r in csv.DictReader(f) if r["source_id"]==SID)
    assert r["source_id"]==SID and len(rows)==4 and len({r["instrument"] for r in rows})==1
    rng=np.random.default_rng(167338)
    t,y,s=np.array([[float(r[k]) for k in ["time_jd","rv_kms","err_kms"]] for r in rows]).T
    t-=DR3_JD
    mu=np.array([float(r[k]) for k in FIELDS])
    L=np.linalg.cholesky(covariance(r))
    models={}
    with (HERE/"desi_lpv_residuals.csv").open("w") as f:
        writer=csv.writer(f)
        writer.writerow(["model","expid","time_jd","rv_kms","err_kms","fitted_rv_kms","residual_kms"])
        for label,scale,extra in [("nominal",1,0.),("inflated",3,5.)]:
            prior=mu+rng.normal(size=(65536,6))@(scale*L.T)
            prior=prior[np.all((prior>GAIA_MIN)&(prior<GAIA_MAX),axis=1)]
            error=np.hypot(s,extra)
            info,idx,fitted=analyze(prior,mu,L*scale,t,y,error,rng,np.ones(len(y),bool))
            info["desi_offset_kms"]=info.pop("lamost_offset_kms")
            info["best"]["desi_offset_kms"]=info["best"].pop("lamost_offset_kms")
            draws=prior[idx]
            fm=1.036149e-7*draws[:,0]*draws[:,2]**3*(1-draws[:,3]**2)**1.5
            info["f_rv_msun"]=np.quantile(fm,[.16,.5,.84]).tolist()
            info["critical_primary_for_m2min_gt3"]=np.quantile(np.sqrt(27/fm)-3,[.16,.5,.84]).tolist()
            m1=float(m["m1"])
            msig=(float(m["m1_upper"])-float(m["m1_lower"]))/2
            info["primary_scenarios"]={}
            for name,mean,sig in [("catalog",m1,msig),("30pct_uncertainty",m1,m1*.3),
                                  ("30pct_lower_primary",m1*.7,0.)]:
                pri=rng.normal(mean,sig,len(fm))
                ok=pri>0
                mm=companion(fm[ok],pri[ok])
                info["primary_scenarios"][name]=dict(mean=mean,sigma=sig,
                         m2min_msun=np.quantile(mm,[.16,.5,.84]).tolist(),p_gt3=np.mean(mm>3))
            models[label]=info
            for row,sig,model,val in zip(rows,error,fitted,y):
                writer.writerow([label,row["expid"],row["time_jd"],val,sig,model,val-model])
    (HERE/"desi_lpv_fits.json").write_text(json.dumps(dict(source_id=SID,span_days=np.ptp(t),
        fraction_period=np.ptp(t)/mu[0],input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [archive,visits,masses]},models=models),indent=2)+"\n")
if __name__=="__main__":
    main()
