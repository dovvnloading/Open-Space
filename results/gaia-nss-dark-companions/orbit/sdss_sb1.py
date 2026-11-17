"""Test SDSS epoch velocities against DR3 SB1 orbits with full Gaia covariance."""
import argparse
import ast
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

from fit_controls import HERE, DR3_JD, TWOPI, companion

FIELDS = ["period","center_of_mass_velocity","semi_amplitude_primary","eccentricity",
          "arg_periastron","t_periastron"]
RAW = HERE/"sdss_sb1_nss.csv"
NAMES = ["1830677037678583680","4661290764764683776","4685960537144695552"]
GAIA_MIN = np.array([.01,-np.inf,.0001,0,-np.inf,-np.inf])
GAIA_MAX = np.array([np.inf,np.inf,np.inf,.999,np.inf,np.inf])


def extract(catalog):
    with gzip.open(catalog,"rt") as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        rows = [r for r in reader if r["source_id"] in NAMES and r["nss_solution_type"]=="SB1"]
    fields = ["source_id","bit_index","corr_vec"]+FIELDS+[k+"_error" for k in FIELDS]
    assert len(rows)==3
    with RAW.open("w") as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def curve(p,t):
    period,gamma,k,e,w,tp = np.atleast_2d(p).T[:,:,None]
    mean = ((TWOPI*(t-tp)/period + np.pi) % TWOPI)-np.pi
    E = mean.copy()
    for _ in range(12):
        step=(E-e*np.sin(E)-mean)/(1-e*np.cos(E))
        E-=step
        if np.max(np.abs(step))<1.e-12:break
    w=np.deg2rad(w)
    return gamma+k*((1-e**2)*np.cos(E)*np.cos(w)-
                    np.sqrt(1-e**2)*np.sin(E)*np.sin(w))/(1-e*np.cos(E))


def covariance(row):
    assert int(row["bit_index"])==127
    C=np.eye(6)
    C[np.tril_indices(6,-1)]=ast.literal_eval(row["corr_vec"])
    C+=np.tril(C,-1).T
    err=np.array([float(row[k+"_error"]) for k in FIELDS])
    return C*np.outer(err,err)


def qtiles(x):
    return np.quantile(x,[.16,.5,.84]).tolist()


def run(sid, row, visits, mass, floor, rng, writer):
    n=len(visits)
    assert len({r["instrument"] for r in visits})==1
    t,y,err=np.array([[float(r[k]) for k in ["time_jd","rv_kms","err_kms"]]
                      for r in visits]).T
    t-=DR3_JD
    s=np.hypot(err,floor)
    mu=np.array([float(row[k]) for k in FIELDS])
    L=np.linalg.cholesky(covariance(row))
    pred=curve(mu,t).ravel()
    prior_draws=mu+rng.normal(size=(16384,6))@L.T
    prior_draws=prior_draws[np.all((prior_draws>GAIA_MIN)&(prior_draws<GAIA_MAX),axis=1)]
    sim=curve(prior_draws,t)+rng.normal(size=(len(prior_draws),n))*s
    # Compare fixed Gaia orbital shape, profiling an arbitrary survey offset.
    # Empirical Monte Carlo rank includes Gaia phase and amplitude errors.
    def shape_stat(velocities):
        dev=velocities-pred
        dev=dev-np.average(dev,axis=-1,weights=1/s**2)[...,None]
        return np.sum((dev/s)**2,axis=-1)
    fixed_chi=float(shape_stat(y))
    p=(1+np.sum(shape_stat(sim)>=fixed_chi))/(len(sim)+1)

    def res(x):
        diff=solve_triangular(L,x[:6]-mu,lower=True,check_finite=False)
        return np.r_[diff, (curve(x[:6],t).ravel()+x[6]-y)/s]
    initial=np.r_[mu,np.average(y-pred,weights=1/s**2)]
    fit=least_squares(res,initial,bounds=(np.r_[GAIA_MIN,-np.inf],np.r_[GAIA_MAX,np.inf]),
                      x_scale="jac",max_nfev=4000)
    assert fit.success
    x=fit.x
    fitted=curve(x[:6],t).ravel()+x[6]
    chi_gaia=float(res(x)[:6]@res(x)[:6])
    chi_rv=float(res(x)[6:]@res(x)[6:])
    info=dict(source_id=sid,used=n,n_quality=sum(int(r["quality"]) for r in visits),
              rv_floor_kms=floor,baseline_days=np.ptp(t),observed_range_kms=np.ptp(y),
              gaia_range_kms=np.ptp(pred),profiled_fixed_chi2=fixed_chi,prior_predictive_p=p,
              prior_predictive_draws=len(sim),joint_chi2_gaia=chi_gaia,joint_chi2_rv=chi_rv,
              joint_dof=n-1,best_fit=dict(zip(FIELDS+["survey_offset_kms"],x)),
              rv_rms_kms=np.sqrt(np.mean((y-fitted)**2)))
    for r,v,pr,ft in zip(visits,y,pred,fitted):
        writer.writerow([sid,r["spectrum_pk"],r["time_jd"],r["quality"],floor,v,
                         r["err_kms"],pr,ft,v-ft])
    # Conditional edge-on mass only. Posterior draws from a local approximation
    # are never an identification of a compact object.
    if mass["m1"]:
        D=np.linalg.cholesky(np.linalg.inv(fit.jac.T@fit.jac))
        z=rng.normal(size=(8192,7))
        pars=x+z@D.T
        good=np.all((pars[:,:6]>GAIA_MIN)&(pars[:,:6]<GAIA_MAX),axis=1)
        pars,z=pars[good],z[good]
        rvdev=(curve(pars[:,:6],t)+pars[:,6,None]-y)/s
        gdev=solve_triangular(L,(pars[:,:6]-mu).T,lower=True).T
        lnw=-.5*((rvdev**2).sum(1)+(gdev**2).sum(1)-(z**2).sum(1))
        weights=np.exp(lnw-lnw.max())
        weights/=weights.sum()
        info["importance_ess"]=1/(weights@weights)
        draws=pars[rng.choice(len(pars),size=8192,p=weights)]
        f=1.036149e-7*draws[:,0]*draws[:,2]**3*(1-draws[:,3]**2)**1.5
        m1=float(mass["m1"])
        sigma=(float(mass["m1_upper"])-float(mass["m1_lower"]))/2
        primaries=rng.normal(m1,sigma,len(f))
        info["conditional_f_rv_msun_q16_q50_q84"]=qtiles(f)
        info["conditional_min_m2_msun_q16_q50_q84"]=qtiles(companion(f[primaries>0],primaries[primaries>0]))
    return info


def main(catalog):
    if catalog:
        extract(catalog)
    eph=HERE.parent/"epochs"/"sdss_visits.csv"
    mpath=HERE.parent/"gaia"/"mass_screen.csv"
    with eph.open() as f: visits=list(csv.DictReader(f))
    with mpath.open() as f:masses={r["source_id"]:r for r in csv.DictReader(f)}
    with RAW.open() as f: gaia={r["source_id"]:r for r in csv.DictReader(f)}
    result=[]
    rng=np.random.default_rng(955104)
    with (HERE/"sdss_sb1_residuals.csv").open("w") as f:
        writer=csv.writer(f)
        writer.writerow(["source_id","spectrum_pk","jd","quality","rv_floor_kms","rv_kms",
                         "formal_error_kms","gaia_prediction_kms","joint_fit_kms","residual_kms"])
        for sid in NAMES:
            vs=[r for r in visits if r["source_id"]==sid]
            # HD 192288 and LHA 120-S 80 are flagged diagnostic fits.
            if sid=="4685960537144695552":
                vs=[r for r in vs if int(r["quality"])]
            for floor in [0.,1.,5.]:
                result.append(run(sid,gaia[sid],vs,masses[sid],floor,rng,writer))
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [RAW,eph,mpath]}
    (HERE/"sdss_sb1_fits.json").write_text(json.dumps(dict(sha256=hashes,fits=result),indent=2)+"\n")

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--full-nss-csv-gz",type=Path)
    main(parser.parse_args().full_nss_csv_gz)
