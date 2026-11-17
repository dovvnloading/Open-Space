"""Joint Keplerian fits to Gaia and published RVs of designated benchmarks."""
import argparse
import ast
import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
DR3_JD = 2457389.0
DR4_JD = 2451545.0 + 17.5 * 365.25
AU_KM = 149597870.7
TWOPI = 2 * np.pi
MASSES = {"BH1": (.93, .05), "BH2": (1.07, .19), "NS1": (.79, .01),
          "BH3": (.76, .05)}
IDS = {"4373465352415301632":"BH1", "5870569352746779008":"BH2",
       "6328149636482597888":"NS1"}


def ti_elements(omega, node, cosi):
    c, s = np.cos(omega), np.sin(omega)
    n, m = np.cos(node), np.sin(node)
    return (n*c-m*s*cosi, m*c+n*s*cosi,
            -n*s-m*c*cosi, -m*s+n*c*cosi)


def ecc_anomaly(time, period, ecc, tp):
    mean = (TWOPI * (time-tp)/period + np.pi) % TWOPI - np.pi
    E = mean.copy()
    for _ in range(12):
        step = (E-ecc*np.sin(E)-mean)/(1-ecc*np.cos(E))
        E -= step
        if np.max(abs(step)) < 1.e-12:
            break
    return E


def rv_model(x, time):
    period, ecc, tp, a, cosi, omega, _, _, gamma = x[:9]
    E = ecc_anomaly(time, period, ecc, tp)
    projected = a * np.sqrt(1-cosi**2)
    speed = TWOPI*projected*AU_KM/(86400*period)
    return gamma + speed * (np.cos(omega)*np.sqrt(1-ecc**2)*np.cos(E)
                            -np.sin(omega)*np.sin(E))/(1-ecc*np.cos(E))


def nss_constraint(row, scale):
    # Gaia DR3 archive data model 20.6.1 lists separate packing orders.
    full = row["nss_solution_type"] == "AstroSpectroSB1"
    if int(row["bit_index"]) != (65535 if full else 8191):
        raise ValueError("Unsupported NSS parameter mask")
    extras = ["c_thiele_innes","h_thiele_innes","center_of_mass_velocity"] if full else []
    fields = (["ra","dec","parallax","pmra","pmdec"] +
              [k+"_thiele_innes" for k in "abfg"] + extras +
              ["eccentricity","period","t_periastron"])
    n = len(fields)
    corr = np.eye(n)
    corr[np.tril_indices(n, k=-1)] = ast.literal_eval(row["corr_vec"])
    corr += np.tril(corr, k=-1).T
    sig = np.array([float(row[k+"_error"]) for k in fields])
    cov = corr * np.outer(sig, sig) * scale**2
    # Extract a marginal block (equivalent to integrating over position/PM).
    keep = [2] + list(range(5,n))
    L = np.linalg.cholesky(cov[np.ix_(keep,keep)])
    obs = np.array([float(row[fields[j]]) for j in keep])
    plx, A, B, F, G = obs[:5]
    u = np.hypot(A+G,B-F)
    v = np.hypot(A-G,B+F)
    alpha = (u+v)/2
    cosi = (u-v)/(u+v)
    ap, am = np.arctan2(B-F,A+G), np.arctan2(B+F,A-G)
    ecc, period, tp = obs[-3:]
    start = np.array([period,ecc,tp,alpha/plx,cosi,(ap-am)/2,(ap+am)/2,plx,
                      float(row["center_of_mass_velocity"] or 0.)])

    def residual(x):
        period,ecc,tp,a,ci,omega,node,parallax,gamma = x[:9]
        abfg = np.array(ti_elements(omega,node,ci))*a*parallax
        extra = []
        if full:
            b = a*np.sqrt(1-ci**2)
            extra = [b*np.sin(omega),b*np.cos(omega),gamma]
        model = np.r_[parallax, abfg, extra, ecc, period, tp]
        return solve_triangular(L,model-obs,lower=True,check_finite=False)
    return residual, start


def black_hole_three(cache):
    dat = np.loadtxt(cache / "epochast.dat")
    dat = dat[dat[:,7] == 0]
    t = dat[:,2]-DR4_JD
    y, error, pf, angle = dat[:,3:7].T
    sina, cosa = np.sin(np.deg2rad(angle)), np.cos(np.deg2rad(angle))
    basis = np.stack([sina,cosa,sina*t/365.25,cosa*t/365.25],axis=1)
    def residual(x):
        period,ecc,tp,a,ci,omega,node,plx = x[:8]
        E = ecc_anomaly(t,period,ecc,tp)
        X, Y = np.cos(E)-ecc, np.sqrt(1-ecc**2)*np.sin(E)
        A,B,F,G = ti_elements(omega,node,ci)
        al = ((B*X+G*Y)*sina+(A*X+F*Y)*cosa)*a*plx
        return (al + pf*plx + basis@x[9:13] - y)/error
    # Approximate published orbit initializes a recovery, not a blind search.
    start = np.array([4200,.72,240,16,-.35,1.36,2.38,1.68,-357,4,2,-28,-155])
    rv = np.loadtxt(cache / "epochrv.dat")
    return residual,start,rv[:,1]-DR4_JD,rv[:,2],rv[:,3]


def companion(f,m1):
    m = f + 2*m1
    for _ in range(12):
        delta = (m**3 - f*(m+m1)**2)/(3*m*m-2*f*(m+m1))
        m -= delta
    return m


def fit(name, start, astro, t, y, error, instrument, rng, output):
    # Parameters use AU, mas, days, radians and km/s. Free nuisance
    # astrometric offsets and proper motions appear only in BH3.
    npar = len(start)
    lower = np.r_[[.1,0,-np.inf,.001,-.99999,-np.inf,-np.inf,.001,-np.inf],
                  np.full(npar-9,-np.inf)]
    upper = np.r_[[1.e5,.99,np.inf,200,.99999,np.inf,np.inf,100,np.inf],
                  np.full(npar-9,np.inf)]
    def resid(x):
        return np.r_[astro(x),(rv_model(x,t)-y)/error]
    fits = []
    for shift in [0,np.pi]:
        x = start.copy()
        x[5:7] += shift
        # Bring precise RVs onto the Gaia phase before restoring formal errors.
        # The period uncertainty can shift the RV predictions by many sigma.
        for inflate in [100., 10., 1.]:
            def tempered(p):
                return np.r_[astro(p), (rv_model(p,t)-y)/(inflate*error)]
            step = least_squares(tempered,x,bounds=(lower,upper),x_scale="jac",
                                 max_nfev=700)
            x = step.x
        fits.append(step)
    result = min(fits,key=lambda f:f.cost)
    if not result.success:
        raise RuntimeError(result.message)
    x = result.x
    residual = y-rv_model(x,t)
    cov = np.linalg.inv(result.jac.T@result.jac)
    # Importance-correct a local multivariate normal proposal. The likelihood
    # uses the actual nonlinear model at every draw, with flat parameter priors.
    L = np.linalg.cholesky(cov)
    z = rng.normal(size=(5000,npar))
    draws = x + z@L.T
    lnw = np.full(len(z),-np.inf)
    for i,p in enumerate(draws):
        if np.all((p>lower)&(p<upper)):
            r = resid(p)
            lnw[i] = -.5*(r@r-2*result.cost-z[i]@z[i])
    w = np.exp(lnw-np.max(lnw))
    w /= w.sum()
    ess = 1/(w@w)
    # Report weighted quantiles. Each Gaussian prior on the luminous mass
    # is independent of the conditional orbital solution.
    m1 = rng.normal(*MASSES[name], len(w))
    valid = (w>0)&(m1>0)
    w, d, m1 = w[valid], draws[valid], m1[valid]
    w /= w.sum()
    f = d[:,3]**3/(d[:,0]/365.25)**2
    sin_i = np.sqrt(1-d[:,4]**2)
    def quant(a):
        order = np.argsort(a)
        weights = w[order]
        return np.interp([.16,.50,.84],np.cumsum(weights)-weights/2,a[order]).tolist()
    fit_info = dict(parameter_order=["period_d","eccentricity","tp_d_rel_ref","a1_AU",
                    "cos_i","omega_rad","Omega_rad","parallax_mas","gamma_kms"]+
                    (["dRA_mas","dDec_mas","pmRA_masyr","pmDec_masyr"] if npar==13 else []),
                    best_fit=x.tolist(), n_rv=len(y),n_gaia_constraints=len(astro(x)),
                    rv_rms_kms=float(np.sqrt(np.mean(residual**2))),
                    rv_chi2=float(np.sum((residual/error)**2)),gaia_chi2=float(astro(x)@astro(x)),
                    joint_dof=len(t)+len(astro(x))-npar, proposal_ess=ess,
                    quantile_levels=[.16,.5,.84],period_d=quant(d[:,0]),ecc=quant(d[:,1]),
                    inclination_deg=quant(np.rad2deg(np.arccos(d[:,4]))),
                    f_rv_msun=quant(f*sin_i**3),m2_min_msun=quant(companion(f*sin_i**3,m1)),
                    m2_msun=quant(companion(f,m1)))
    for time, inst, v, s, dev in zip(t,instrument,y,error,residual):
        output.writerow([name,f"{time:.6f}",inst,v,s,f"{v-dev:.6f}",f"{dev:.6f}"])
    return fit_info


def main(cache):
    rng = np.random.default_rng(110572)
    with (HERE/"nss_controls.csv").open() as f:
        gaia = {IDS[r["source_id"]]:r for r in csv.DictReader(f)}
    with (HERE/"published_rv.csv").open() as f:
        epochs = list(csv.DictReader(f))
    info = {}
    with (HERE/"control_residuals.csv").open("w") as f:
        writer = csv.writer(f)
        writer.writerow(["control","jd_rel_reference","instrument","rv_kms","err_kms",
                         "model_kms","observed_minus_model_kms"])
        for name in MASSES:
            if name == "BH3":
                astro,start,t,y,s = black_hole_three(cache)
                inst = ["Gaia"]*len(t)
            else:
                rows = [r for r in epochs if r["control"]==name]
                t,y,s = np.array([[float(r[k]) for k in ["hjd_utc","rv_kms","err_kms"]]
                                  for r in rows]).T
                t -= DR3_JD
                inst = [r["instrument"] for r in rows]
                astro,start = nss_constraint(gaia[name],3. if name=="NS1" else 1.)
                if name!="BH2":
                    start[8] = np.median(y)
            info[name] = fit(name,start,astro,t,y,s,inst,rng,writer)
            print(name,info[name]["m2_msun"],info[name]["proposal_ess"], flush=True)
    (HERE/"control_fits.json").write_text(json.dumps(info,indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cache",type=Path,help="Directory with ESA BH3 timeseries")
    main(parser.parse_args().cache)
