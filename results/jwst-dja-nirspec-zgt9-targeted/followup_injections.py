#!/usr/bin/env python3
"""Reproducible *limited* per-target interloper injection / recovery diagnostic.

Tests the exact PR #195 1D fitter on independently parametrized dusty-Balmer,
line-confusion, cool-star and artefact surrogate spectra using target-specific
wave/error/valid masks. Samples are NOT an empirical astrophysical prior or
calibrated redshift posterior; don't promote conditional high-z quantiles.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from astropy.io import fits

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import run_dispatch as rd
import download_targets as down

WORK=Path("/tmp/dja-pr197-nine-fits")
LEDGER=HERE/"literature_referee_2026-11-19.json"

def json_ready(x):
    if isinstance(x,dict):return {str(k):json_ready(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [json_ready(v) for v in x]
    if isinstance(x,np.ndarray):return json_ready(x.tolist())
    if isinstance(x,np.generic):return json_ready(x.item())
    if isinstance(x,float) and not np.isfinite(x):return None
    return x

def line(w,rest,z,R):
    centre=rest*(1+z)
    sigma=max(centre/R/2.355,.0008)
    return np.exp(-.5*((w-centre)/sigma)**2)

def template(w,family,z,rng,band):
    # This generator is separate from and deliberately broader than the fitted
    # finite template bank. It is NOT a physically exhaustive SED simulator.
    if family=="dusty_balmer":
        br=.400*(1+z)
        blue=rng.uniform(.015,.55)
        dust=rng.uniform(.25,2.2)
        slope=rng.uniform(-.8,2.)
        continuum=(w/2.)**slope*np.exp(-dust*(2./np.maximum(w,.2))**.7)
        return continuum*np.where(w<br,blue,1.)+rng.uniform(0,1.4)*line(w,.5007,z,130 if band=="prism" else 1000)
    if family=="line_confusion":
        slope=rng.uniform(-1.,1.)
        cont=(w/2.)**slope * rng.uniform(.01,.25)
        k=rng.choice([.3727,.4861,.5007,.6563])
        return cont + rng.uniform(2,12)*line(w,k,z,100 if band=="prism" else 1000)
    if family=="cool_dwarf_surrogate":
        slope=rng.uniform(-2,1.7)
        v=(w/2.)**slope
        for c in [1.15,1.40,1.90,2.70,3.30,4.30]:
            width=rng.uniform(.035,.20)
            depth=rng.uniform(.2,.85)
            v*=1-depth*np.exp(-.5*((w-c)/width)**2)
        return v
    if family=="detector_spike_surrogate":
        cont=np.ones_like(w)*rng.uniform(.0,.3)
        idx=int(rng.integers(3,len(w)-3))
        cont[idx:idx+2]+=rng.uniform(4,12)
        return cont
    if family=="highz_recovery":
        lya=.121567*(1+z)
        return (w/2)**rng.uniform(-1,1)*np.where(w<lya,rng.uniform(0,.025),1.)
    raise ValueError(family)

def one(g,s,model,n_per_family):
    root,name,url=down.safe(s)
    path=WORK/root/name
    if not path.exists():
        down.download_one({"root":root,"file":name,"url":url},WORK,cap_mb=400)
    with fits.open(path) as h:
        w,f,e,valid=rd.arrays_1d(h)
    band="prism" if "prism-clear" in name else "grating"
    bank=model.bank(w,band)
    seed=20261119+int(g["group_id"])*101+sum(map(ord,name))
    rng=np.random.default_rng(seed)
    valid=valid & np.isfinite(f)&np.isfinite(e)&(e>0)
    wt=np.where(valid,1/np.maximum(e,1e-10)**2,0)
    noise_err=np.where(valid,e,1.)
    s2=float(np.sum((np.where(valid,f,0)/noise_err)**2))
    approx_net_sn=float(np.sqrt(max(0,s2-int(valid.sum()))))
    target_sn=float(np.clip(approx_net_sn,4,35))
    families=["dusty_balmer","line_confusion","cool_dwarf_surrogate",
              "detector_spike_surrogate","highz_recovery"]
    spectra=[]; meta=[]
    for family in families:
        for j in range(n_per_family):
            if family=="dusty_balmer": z=rng.uniform(1.5,5.6)
            elif family=="line_confusion": z=rng.uniform(.5,8.)
            elif family=="highz_recovery": z=rng.uniform(12.1,17.)
            else: z=float("nan")
            t=template(w,family,z,rng,band)
            shape_sn=float(np.sqrt(np.sum(t*t*wt)))
            amp=target_sn*rng.uniform(.55,1.45)/max(shape_sn,1e-10)
            mock=amp*t+rng.normal(0,noise_err)
            spectra.append(np.where(valid,mock,0))
            meta.append((family,z))
    # Include original target: all columns pass through exactly the same
    # bank(w, disperser) + fit(w,f,e,valid,bank) as the discovery proxy.
    spectra.insert(0,np.where(valid,f,0.))
    n=len(spectra)
    flux=np.column_stack(spectra)
    errors=np.broadcast_to(noise_err[:,None],(len(w),n)).copy()
    v=np.broadcast_to(valid[:,None],(len(w),n)).copy()
    fit=model.fit(w,flux,errors,v,bank)
    classes={}
    for family in families:
        cols=[j+1 for j,(k,z) in enumerate(meta) if k==family]
        high=np.asarray(fit["p_gt12_proxy"])[cols]
        ok=np.asarray(fit["two_sided"])[cols]
        count=int(np.sum((high>=.5)&ok))
        classes[family]={"n":len(cols),"false_highz_threshold_and_two_sides":count,
            "rate":count/len(cols),"proxy_median":float(np.median(high)),
            "proxy_max":float(np.max(high))}
    result={
        "group_id":g["group_id"],"file":name,"sha256":down.sha256(path),
        "seed":seed,"disperser_fit":band,"valid_bins":int(valid.sum()),
        "effective_total_sn_estimate":approx_net_sn,"injection_target_snr":target_sn,
        "sample_size_per_family":n_per_family,
        "observed_proxy_p_gt12":float(fit["p_gt12_proxy"][0]),
        "observed_p_gt14_proxy":float(fit["p_gt14_proxy"][0]),
        "observed_two_sided":bool(fit["two_sided"][0]),
        "conditional_highz_z16":float(fit["q16_high"][0]),
        "conditional_highz_z50":float(fit["q50_high"][0]),
        "conditional_highz_z84":float(fit["q84_high"][0]),
        "injection_recovery":classes,
        "physical_posterior":None,
        "interpretation":"Synthetic diagnostic FPR on same 1D fitting pipeline; NOT a calibrated physical redshift posterior. No native DQ or NIRCam.",
    }
    return result

def run(out,n):
    out.mkdir(parents=True,exist_ok=True)
    ledger=json.loads(LEDGER.read_text())
    groups=[g for g in ledger["entries"] if g["referee_class"]=="UNRESOLVED_NOT_NEW"]
    model=rd.load_model()
    r={"method":"same PR195 1D proxy, independent surrogate 1D injections, fixed RNG",
       "status":"DIAGNOSTIC_INJECTIONS_NOT_CALIBRATED_POSTERIORS",
       "n_injections_per_family_per_spectrum":n,
       "spectra":[]}
    for g in groups:
        for s in g["individual_spectra"]:
            try:
                result=one(g,s,model,n)
                r["spectra"].append(result)
            except Exception as exc:
                r["spectra"].append({"group_id":g["group_id"],"file":s["file"],
                    "status":"INJECTION_FAILURE", "error":type(exc).__name__+": "+str(exc)[:350]})
            (out/"followup_injections.json").write_text(json.dumps(json_ready(r),indent=2,allow_nan=False)+"\n")
    r["group_count"]=len(groups)
    r["spectra_completed"]=sum("injection_recovery" in s for s in r["spectra"])
    r["any_physical_redshift_posteriors"]=False
    (out/"followup_injections.json").write_text(json.dumps(json_ready(r),indent=2,allow_nan=False)+"\n")
    print("INJECTION_DIAGNOSTIC "+json.dumps({"groups":r["group_count"],"spectra":r["spectra_completed"]}))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--per-family",type=int,default=20)
    a=p.parse_args()
    if not 5<=a.per_family<=100: p.error("per-family must be 5-100")
    run(a.out,a.per_family)
