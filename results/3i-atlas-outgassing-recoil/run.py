#!/usr/bin/env python3
"""Reproduce a new 3I full-arc astrometric fit from live MPC on a networked runner.

Uses the repository's pre-existing independently written N-body sky-plane fitter;
does NOT repeat its incoming-hyperbola classification. All observatory/weighting
and outlier selections are frozen in README before the present run.
"""
from __future__ import annotations
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import inspect
import io
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
import traceback
import urllib.request
import urllib.error

import numpy as np
from jplephem.spk import SPK

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "generated"
OUT.mkdir(exist_ok=True)
OLD = ROOT / "results/inbound-hyperbolic-search/orbits"
sys.path.insert(0, str(OLD))
import acquire
import refit

SCRATCH = Path(tempfile.mkdtemp(prefix="3i-recoil-"))
acquire.CACHE = SCRATCH
refit.CACHE = SCRATCH
refit.DATA = SCRATCH
receipts = []
now = lambda: datetime.now(timezone.utc).isoformat()

def sha(data):
    return hashlib.sha256(data).hexdigest()

def write(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False, default=float) + "\n")

def fetch(url, label, body=None):
    req = urllib.request.Request(url, data=body, method="GET", headers={
        "User-Agent": "Open-Space/3I-recoil-2026 (public scientific research)",
        **({"Content-Type": "application/json"} if body else {}),
    })
    item = {"source":label, "url":url, "request_utc":now()}
    receipts.append(item)
    try:
        with urllib.request.urlopen(req, timeout=110) as response:
            raw = response.read()
            item.update(status=response.status, final_url=response.url, bytes=len(raw),
                        sha256=sha(raw), last_modified=response.headers.get("Last-Modified"))
        return raw
    except Exception as e:
        item["error"] = str(e)
        raise

def mpc_rows():
    # Frozen previous independent orbit analysis supplies an *auditable fallback*,
    # not a pretense of fresh completeness.
    frozen = gzip.decompress((OLD / "data/3I.obs.gz").read_bytes()).decode("ascii", "replace").splitlines()
    candidates = [("prior_MPC_frozen_2026_11_17", frozen)]
    for url in [
        "https://www.minorplanetcenter.net/iau/ECS/MPCAT-OBS/CmtObs.txt.gz",
        "https://minorplanetcenter.net/iau/ECS/MPCAT-OBS/CmtObs.txt.gz",
    ]:
        try:
            data = fetch(url, "MPC comet 80-column bulk")
            lines = acquire.select(gzip.decompress(data).decode("ascii", "replace").splitlines(),
                                   ["0003I", "K25N010"])
            if lines:
                candidates.append(("live_MPC_bulk", lines))
                break
        except Exception as e:
            print("MPC bulk failure:", str(e), flush=True)
    # Ask the live object service as a separate completeness check. MPC currently
    # accepts the documented JSON body on GET for output_format=OBS80.
    try:
        data = fetch("https://data.minorplanetcenter.net/api/get-obs", "MPC get-obs 3I",
                     json.dumps({"desigs":["3I"],"output_format":["OBS80"]}).encode())
        response = json.loads(data)
        rows = response[0]["OBS80"]
        lines = acquire.select(rows.splitlines() if isinstance(rows,str) else rows,
                               ["0003I","K25N010"])
        if lines:
            candidates.append(("live_MPC_get_obs", lines))
    except Exception as e:
        print("MPC get-obs failure:", str(e), flush=True)
    scores = [(len(set(lines)), max((l[15:32] for l in lines),default=""), label) for label,lines in candidates]
    # First optimize latest observation, then distinct-row count; retain all rows
    # within the chosen source, never silently concatenate distinct reductions.
    scored = sorted(candidates, key=lambda q:(max((s[15:32] for s in q[1]),default=""),len(set(q[1]))))
    label,lines = scored[-1]
    data = ("\n".join(lines)+"\n").encode("ascii", "replace")
    (SCRATCH / "3I.obs").write_bytes(data)
    (OUT / "mpc_3i.obs.gz").write_bytes(gzip.compress(data,mtime=0))
    return {"chosen":label,"candidate_counts":scores,"raw_80_column_sha256":sha(data),
            "rows":len(lines),"first":min(s[15:32] for s in lines),
            "last":max(s[15:32] for s in lines),
            "frozen_baseline_rows":len(frozen),
            "scope_note":"MPC ground optical filters and clipping apply downstream; all submitted modes remain frozen"}

def prepare():
    input_info = mpc_rows()
    data = fetch("https://ssd-api.jpl.nasa.gov/sbdb.api?spk=1004083&full-prec=true",
                 "JPL SBDB published fit")
    sbdb=json.loads(data)
    write("sbdb.json", sbdb)
    input_info["jpl_orbit_id"] = sbdb.get("orbit",{}).get("orbit_id")
    input_info["jpl_last_obs"] = sbdb.get("orbit",{}).get("last_obs")
    input_info["jpl_first_obs"] = sbdb.get("orbit",{}).get("first_obs")
    input_info["jpl_n_obs_used"] = sbdb.get("orbit",{}).get("n_obs_used")
    obs = None
    for url in ["https://minorplanetcenter.net/iau/lists/ObsCodes.html",
                "https://www.projectpluto.com/obscodes.htm"]:
        try:
            obs=fetch(url,"MPC observatory positions")
            if b"Greenwich" not in obs: raise ValueError("no Greenwich station")
            (SCRATCH/"stations.txt").write_bytes(obs)
            break
        except Exception as e:
            print("station URL failure:",e,flush=True)
    if not (SCRATCH/"stations.txt").exists():
        raise RuntimeError("Cannot retrieve MPC observatory positions")
    kernel=fetch("https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp",
                 "DE440s planetary kernel")
    expected="c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2"
    if sha(kernel)!=expected: raise RuntimeError("DE440s SHA mismatch")
    (SCRATCH/"de440s.bsp").write_bytes(kernel)
    write("inputs.json", input_info)
    return sbdb

def inverse_square_derivative():
    """Substitute g(r)=r^-2, dg/dr=-2r^-3 into the frozen analytic RTN variational equations."""
    src=inspect.getsource(refit.derivative)
    original='''        z=(r/2.808)**5.093
        g=.111262*(r/2.808)**(-2.15)*(1+z)**(-4.6142)
        dg_dr=g*(-2.15-4.6142*5.093*z/(1+z))/r'''
    replacement='''        g=1./r**2
        dg_dr=-2.*g/r'''
    if src.count(original)!=1: raise RuntimeError("Upstream Marsden code changed; refusing unreviewed substitution")
    environment={}
    exec(src.replace(original,replacement),refit.__dict__,environment)
    return environment["derivative"]

def residual(x, epoch, ts, observers, east, north, positions):
    path=refit.integrate(x,epoch,min(ts)-refit.OBS_TIME_LIMIT,positions).sol
    emit=ts.copy()
    for k in range(5):
        delta=path(emit)[:3].T-observers
        distance=np.linalg.norm(delta,axis=1)
        emit=ts-distance/refit.C
    delta=path(emit)[:3].T-observers
    u=delta/np.linalg.norm(delta,axis=1)[:,None]
    return refit.ARCSEC * np.column_stack(
       [np.einsum("ij,ij->i",u,east),np.einsum("ij,ij->i",u,north)])

def model_fit(name, x0, epoch, t, obs, east, north, sigma, positions, fitmask,
              derivative, inherited):
    refit.derivative=derivative
    n=len(x0)
    state,cov,mask,rr,redchi=refit.fit(x0,epoch,t[fitmask],obs[fitmask],
        east[fitmask],north[fitmask],sigma[fitmask],positions)
    all_res=residual(state,epoch,t,obs,east,north,positions)
    good=np.flatnonzero(fitmask)[mask]
    omitted=np.setdiff1d(np.flatnonzero(fitmask),good)
    result={"force_law":name,"epoch_jd_tdb":epoch,"state_ssb_icrf":state.tolist(),
            "fit_count":len(good),"excluded_by_clipping":omitted.tolist(),
            "weighted_reduced_chi2":float(redchi),
            "rms_fit_per_coordinate_arcsec":float(np.sqrt(np.mean(all_res[good]**2))),
            "rms_all_per_coordinate_arcsec":float(np.sqrt(np.mean(all_res**2))),
            "covariance":cov.tolist(),
            "weight_rule":"1 arcsec per coordinate; station-night effective max 4; 4-sigma two-coordinate rejection",
            "radial_transverse_normal_A_1e8_au_d2":state[6:9].tolist() if n==9 else None}
    if n==9:
        result["A_formal_1sigma_1e8_au_d2"]=np.sqrt(np.diag(cov)[6:9]).tolist()
        result["A_3d_norm_1e8_au_d2"]=float(np.linalg.norm(state[6:9]))
    if not fitmask.all():
        keptout=np.flatnonzero(~fitmask)
        result["withheld_count"]=len(keptout)
        result["withheld_rms_per_coord_arcsec"]=float(np.sqrt(np.mean(all_res[keptout]**2)))
    return result

def run_orbit():
    old_water=refit.derivative
    r2=inverse_square_derivative()
    fits={}
    with SPK.open(str(SCRATCH/"de440s.bsp")) as eph:
        t,observers,east,north,sigma,groups,dropped=refit.observations("3I",eph)
        epoch=float(np.ceil(t.max())+1)
        positions=refit.ephemeris(eph,float(t.min())-refit.OBS_TIME_LIMIT-2,
                                  epoch+1,.5)
        initial=refit.seed("3I",epoch,receipts)
        summary={"n_supported":len(t),"n_unique_station_night":len(set(groups)),
                 "jd_tdb_min":float(t.min()),"jd_tdb_max":float(t.max()),
                 "excluded_non_optical_or_unlocated":dropped,"fit_epoch":epoch,
                 "initial_jpl_seed_state":initial.tolist(),
                 "observed_before_perihelion":int((t<2460977.5).sum()),
                 "observed_after_perihelion":int((t>=2460977.5).sum()),
                 "frozen_source_file_sha256":sha((SCRATCH/"3I.obs").read_bytes())}
        # October 29 2025 00 UTC is JD 2460977.5; perihelion date is not
        # used for selection, weights or orbit initialization.
        valid=np.ones(len(t),dtype=bool)
        seeds=initial.copy()
        for law in ["gravity","r2_rtn","marsden_water_rtn"]:
            try:
                derivative=old_water if law!="r2_rtn" else r2
                p=np.r_[seeds,np.zeros(3)] if law!="gravity" else seeds
                fits[law]=model_fit(law,p,epoch,t,observers,east,north,
                                    sigma,positions,valid,derivative,old_water)
                if law=="gravity": seeds=np.array(fits[law]["state_ssb_icrf"])
                print("FIT",law,fits[law]["fit_count"],
                      fits[law]["rms_fit_per_coordinate_arcsec"],
                      fits[law]["radial_transverse_normal_A_1e8_au_d2"],flush=True)
            except Exception as e:
                fits[law]={"error":str(e),"traceback":traceback.format_exc()[-2500:]}
                print("FIT FAILED",law,str(e),flush=True)
        # Temporal validation is genuinely withheld: latest ten percent by
        # receipt epoch, *not* included in fitting or clipping.
        if "state_ssb_icrf" in fits.get("r2_rtn",{}):
            try:
                cut=np.quantile(t,.90)
                fitmask=t<cut
                hold=model_fit("r2_rtn_last10percent_temporal_holdout",
                    np.array(fits["r2_rtn"]["state_ssb_icrf"]),epoch,t,
                    observers,east,north,sigma,positions,fitmask,r2,old_water)
                fits["temporal_holdout"]=hold
                print("HOLDOUT",hold["withheld_count"],
                      hold["withheld_rms_per_coord_arcsec"],flush=True)
            except Exception as e:
                fits["temporal_holdout"]={"error":str(e)}
        summary["fit_models"]=fits
        write("orbit.json",summary)
    if "state_ssb_icrf" not in fits.get("r2_rtn",{}):
        raise RuntimeError("Independent three-direction inverse-square orbit fit did not succeed")

def main():
    try:
        prepare()
        run_orbit()
        import subprocess
        subprocess.run([sys.executable,str(HERE/"recoil.py")],check=True)
    finally:
        write("retrieval_receipts.json",receipts)
        write("pipeline.json",{"retrieved_utc":now(),"run_py_sha256":sha(Path(__file__).read_bytes()),
            "recoil_py_sha256":sha((HERE/"recoil.py").read_bytes()),
            "production_csv_sha256":sha((HERE/"production.csv").read_bytes()),
            "prior_orbit_refit_py_sha256":sha((OLD/"refit.py").read_bytes()),
            "kernel_expected_sha256":"c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2"})
if __name__=="__main__":
    main()
