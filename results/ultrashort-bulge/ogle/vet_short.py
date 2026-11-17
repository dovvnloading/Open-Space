"""Local multi-start fits and baseline diagnostics for all four EWS sub-half-day BLG alerts.

Run python vet_short.py CACHE. VBMicrolensing, NumPy, SciPy, and pyerfa required.
"""
import argparse
import csv
import json
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import VBMicrolensing
from scipy.optimize import differential_evolution, minimize
from scipy.stats import spearmanr

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "ogle-dark-lenses"))
import lens
from freeze import catalog, fetch_curve


def boxfit(fn, bounds, seed, iterations, starts=()):
    opt = differential_evolution(fn, bounds, seed=seed, maxiter=iterations, popsize=7,
                                 polish=False, tol=0.0001)
    best = opt.fun, opt.x
    for p in [opt.x, *starts]:
        # EWS can print a zero impact parameter for high magnification.
        p=np.clip(p,[lo for lo,_ in bounds],[hi for _,hi in bounds])
        result = minimize(fn, p, method="Nelder-Mead", bounds=bounds,
                          options=dict(maxiter=1100, xatol=0.0001, fatol=0.005))
        if result.fun < best[0]:
            best = result.fun, result.x
    return best


def solve(arguments):
    cache, row = arguments
    y, _, n = row["event"].split("-")
    path = cache / y / f"blg-{n}.dat"
    t, mag, err, see, sky = lens.read_phot(path)
    t -= float(row["t0_hjd"])-2450000
    flux, sig = lens.to_flux(mag, err)
    valid = np.isfinite(sig) & (sig>0)
    off = (np.abs(t)>30) & valid
    base = np.median(flux[off])
    s = max(1.0, np.median(np.abs((flux[off]-base)/sig[off]))/0.67449)
    rn = (flux-base)/(sig*s)
    outlier_groups = []
    civil_night = np.floor(t + float(row["t0_hjd"]))  # noon UT boundary
    for night in np.unique(civil_night[off]):
        mask = off & (civil_night==night)
        if np.sum(mask & (rn>5))>=3:
            outlier_groups.append(dict(hjd_noon=int(night), n=int(mask.sum()), median_z=float(np.median(rn[mask]))))
    near = (np.abs(t)<15) & valid
    times = t[near].tolist()
    f, sig = flux[near], sig[near]*s
    w = 1/sig**2
    base_w = (w*f).sum()/w.sum()
    def fitflux(amp, details=False):
        # Weighted least-squares source and blend, with each flux nonnegative.
        amp = np.asarray(amp)
        fs, fb = lens.linear_fluxes(amp, f, sig)
        if fs <= 0:
            fs, fb = 0., base_w
        elif fb < 0:
            fs, fb = max(0.,(w*amp*f).sum()/(w*amp**2).sum()), 0.
        pred = amp*fs+fb
        stat = np.sum(((f-pred)/sig)**2)
        if not np.isfinite(stat):
            stat = 1e50
        return (float(stat), float(fs), float(fb)) if details else stat

    vb = VBMicrolensing.VBMicrolensing()
    vb.LoadESPLTable(str(pathlib.Path(VBMicrolensing.__file__).parent / "data/ESPL.tbl"))
    vb.Tol = 0.0001
    def amp(p, model):
        dt, lu, lt = p[:3]
        if model=="flare_pulse":
            offsets=t[near]-dt
            # Independent exponential rise/decay e-folding times with a fast-rise constraint.
            rise, decay=np.exp(lu),np.exp(lu)+np.exp(lt)
            return 1+10*np.exp(np.where(offsets<0,offsets/rise,-offsets/decay))
        if model=="PSPL":
            return lens.magnification(np.exp(2*lu)+((t[near]-dt)/np.exp(lt))**2)
        if model=="FSPL":
            return vb.ESPLLightCurve([lu,lt,dt,p[3]],times)[0]
        if model=="asymmetric":
            te = np.where(t[near]<dt, np.exp(lt),np.exp(p[3]))
            return lens.magnification(np.exp(2*lu)+((t[near]-dt)/te)**2)
        return vb.BinaryLightCurve([p[3],p[4],np.exp(lu),p[5],p[6],lt,dt],times)[0]
    def logbounds(lo,hi):
        return (np.log(lo),np.log(hi))

    common = [(-1.5,1.5),logbounds(.001,5),logbounds(.01,150)]
    allbounds = dict(PSPL=common, FSPL=common+[logbounds(.003,15)],
                     asymmetric=common+[logbounds(.003,150)],
                     flare_pulse=[(-1.5,1.5),logbounds(.001,2),logbounds(.0001,5)],
                     binary=common+[logbounds(.35,3.5),logbounds(.0003,1),
                                    (0,2*np.pi),logbounds(.0003,.2)])
    np.seterr(all="ignore")
    fits = {}
    for kind, bounds in allbounds.items():
        fn=lambda p:fitflux(amp(p,kind))
        starts=[]
        if kind=="PSPL":
            starts=[[0,np.log(max(.001,float(row["umin"]))),np.log(float(row["tau_days"]))]]
        else:
            single=fits["PSPL"]["params"]
            if kind=="FSPL":
                starts=[single+[np.log(r)] for r in [.01,.3,1,3]]
            if kind=="asymmetric":
                starts=[single+[single[2]]]
            if kind=="binary":
                starts=[single+[np.log(s0),np.log(q),ang,np.log(.01)]
                        for s0,q,ang in [(1.2,.01,1.7),(.7,.3,4),(2.1,.1,2),(1.3,.001,5)]]
        stat, p = boxfit(fn,bounds,301,75 if kind=="binary" else 50,starts)
        chi,fs,fb=fitflux(amp(p,kind),True)
        fits[kind]=dict(chi2=chi, params=p.tolist(), t0_offset=float(p[0]),
                         u0=float(np.exp(p[1])), tE_days=float(np.exp(p[2])), fs=fs, fb=fb)
        if kind=="flare_pulse":
            fits[kind].pop("u0"); fits[kind].pop("tE_days")
            fits[kind].update(rise_e_days=float(np.exp(p[1])),decay_e_days=float(np.exp(p[1])+np.exp(p[2])))
        if kind=="FSPL":
            fits[kind]["rho"]=float(np.exp(p[3]))
        if kind=="binary":
            fits[kind].update(s=float(np.exp(p[3])),q=float(np.exp(p[4])),rho=float(np.exp(p[6])))
        if kind=="asymmetric":
            fits[kind]["fall_tE_days"]=float(np.exp(p[3]))
        if kind in ("PSPL","FSPL"):
            restrictions={"long":(2,logbounds(.5,150)),"short":(2,logbounds(.01,.5))}
            if kind=="FSPL":
                restrictions.update(rho_gt1=(3,logbounds(1.,15)),rho_lt1=(3,logbounds(.003,1.)))
            for label,(idx,lim) in restrictions.items():
                if lim[0]<=p[idx]<=lim[1]:
                    fits[kind][label+"_chi2"]=chi
                    continue
                restricted=list(bounds);restricted[idx]=lim
                x=p.copy();x[idx]=np.clip(x[idx],*lim)
                trial,_=boxfit(fn,restricted,709,45,[x])
                fits[kind][label+"_chi2"]=float(trial)
    a = np.sort(t[(abs(t)<=1) & valid])
    pre = a[a<=0];post = a[a>0]
    res=dict(event=row["event"], assumed_t0_hjd=float(row["t0_hjd"]), n_curve=len(t), n_local=int(near.sum()),
             n_before_30d=int(np.sum(t< -30)),n_after_30d=int(np.sum(t>30)), baseline_years=round(float(np.ptp(t)/365.25),2),
             n_halfday=int(np.sum(abs(a)<=.5)), last_pre_peak_offset_days=float(pre[-1]) if len(pre) else None,
             first_post_peak_offset_days=float(post[0]) if len(post) else None, error_scale=float(s),
             baseline_median_i=float(np.median(mag[off])),
             baseline_r_spearman_seeing=float(spearmanr(mag[off],see[off]).statistic),
             baseline_r_spearman_sky=float(spearmanr(mag[off],sky[off]).statistic),
             baseline_nights_at_least_three_5sigma_high=outlier_groups,
             chi2_constant=fitflux(np.ones_like(f)),fits=fits)
    return res


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=pathlib.Path)
    parser.add_argument("--events", nargs="+", help="Specific additional BLG IDs to fit")
    parser.add_argument("--output", default="short-fits.json")
    args=parser.parse_args()
    rows=list(csv.DictReader((HERE/"catalog.csv").open()))
    if args.events:
        other_years={e.split('-')[0] for e in args.events}-{"2025","2026"}
        for yr in sorted(other_years):
            _,season=catalog(args.cache,yr)
            for row in season:
                if row["event"] in args.events:
                    state=fetch_curve(args.cache,row)
                    if state["status"]!="ok":raise RuntimeError(state)
            rows.extend(season)
        short=[r for r in rows if r["event"] in args.events]
        if len(short)!=len(set(args.events)): raise ValueError("unknown event ID")
    else:
        short=[r for r in rows if 0<float(r["tau_days"])<.5]
    with ProcessPoolExecutor(4) as pool:
        records=list(pool.map(solve,[(args.cache,r) for r in short]))
    (HERE/args.output).write_text(json.dumps(records,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
