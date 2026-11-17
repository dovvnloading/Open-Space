"""Independent 250 au refits of the fourteen timed-out, station-night-sampled recent arcs."""
import base64
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
import sys
import traceback

import numpy as np
from jplephem.spk import SPK

ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent / "slow"))
import replay as slow
from replay import acquire, fit

DATA=ROOT / "evidence"
CACHE=ROOT / "cache"
DATA.mkdir(exist_ok=True)
CACHE.mkdir(exist_ok=True)
acquire.CACHE=fit.CACHE=CACHE
fit.DATA=DATA
INPUT=ROOT.parents[1] / "batch/evidence/recent_stratified.json.gz"

def fit_record(row):
    spk=row["spkid"]
    name=row["request"]["desigs"][0]
    obj=name.split("/")[-1].replace(" ","-")
    record={"name":row["name"],"spkid":spk,"queries":[],"models":{}}
    try:
        raw=("\n".join(row["selected_mpc80"])+"\n").encode("ascii")
        sha=hashlib.sha256(raw).hexdigest()
        if sha!=row["fit"]["astrometry_sha256"]:
            raise ValueError("Frozen selection hash mismatch")
        (DATA/(obj+".obs")).write_bytes(raw)
        record.update(selected=len(row["selected_mpc80"]),selected_sha256=sha)
        with SPK.open(str(CACHE/"de440.bsp")) as eph:
            t,observers,east,north,sigma,groups,dropped=fit.observations(obj,eph)
            epoch=float(np.ceil(t.max())+1)
            initial=slow.initial_state(obj,spk,epoch,eph,record["queries"])
            positions=fit.ephemeris(eph,t.min()-fit.OBS_TIME_LIMIT-2,epoch+1,.25)
            record.update(epoch_jd_tdb=epoch,initial_state_ssb_icrf=initial.tolist(),
                          unique_station_nights=len(set(groups)),supported=len(t),
                          supported_span_days=float(t.max()-t.min()),dropped=dropped)
            seed=initial
            modes=["gravity"]
            for mode in modes:
                res={}
                record["models"][mode]=res
                try:
                    ng={"gravity":0,"radial_r2":1,"marsden_outward":2}[mode]
                    x,cov,mask,r,redchi,predict=fit.fit(np.r_[seed,np.zeros(ng)],epoch,t,observers,
                                          east,north,sigma,positions,outward=mode=="marsden_outward",
                                          return_predictor=True)
                    if mode=="gravity": seed=x.copy()
                    res.update(used=int(mask.sum()),rms_per_coordinate_as=float(np.sqrt(np.mean(r[mask]**2))),
                               weighted_reduced_chi2=float(redchi),state_ssb_icrf=x.tolist(),
                               covariance=cov.tolist())
                    res["boundaries"]=fit.inbound(obj,x,cov,epoch,eph)
                    b=res["boundaries"].get("250au") or {}
                    alpha=b.get("inverse_a_au_inv",0.)
                    s=b.get("formal_sigma_inverse_a",np.inf)
                    if mode=="gravity" and alpha < -2.5*s and spk!="1004083":
                        modes.extend(["radial_r2","marsden_outward"])
                        res["holdouts"]=fit.check_holdouts(obj,x,epoch,eph,t,observers,east,north,sigma,
                                                         groups,positions)
                    if alpha<0 and mode!="gravity":
                        try: res["zero_energy_profile"]=fit.parabolic_profile(
                                  epoch,eph,x,mode,mask,sigma,predict)
                        except Exception as exc: res["zero_energy_profile"]={"failed":str(exc)}
                except Exception as exc:
                    res["failed"]=str(exc)
                    if mode=="gravity" and spk!="1004083":
                        modes.extend([m for m in ("radial_r2","marsden_outward") if m not in modes])
    except Exception as exc:
        record["error"]=str(exc)
        traceback.print_exc()
    print("RECENT_RESULT",spk,json.dumps(record,allow_nan=False),flush=True)
    return record

def main():
    source=INPUT.read_bytes()
    rows=[r for r in json.loads(gzip.decompress(source))
          if r.get("fit",{}).get("status")=="fit timeout" and r.get("selected_mpc80")]
    if len(rows)!=14: raise ValueError("Expected fourteen frozen timed-out selections")
    receipts=[]
    kernel=acquire.retrieve(slow.PLANETS,"de440.bsp",receipts)
    if hashlib.sha256(kernel).hexdigest()!=slow.KERNEL_SHA: raise ValueError("Invalid DE440 checksum")
    if b"Greenwich" not in acquire.retrieve(slow.STATIONS,"stations.txt",receipts):
        raise ValueError("Bad observatory table")
    with ProcessPoolExecutor(max_workers=3) as pool: result=list(pool.map(fit_record,rows))
    out={"frozen_archive_sha256":hashlib.sha256(source).hexdigest(),"inputs":receipts,
         "objects":{r["spkid"]:r for r in result}}
    b=(json.dumps(out,indent=2,allow_nan=False)+"\n").encode()
    (DATA/"results.json").write_bytes(b)
    compressed=base64.b64encode(gzip.compress(b,mtime=0)).decode()
    for i in range(0,len(compressed),6000):
        print("RECENT_ARCHIVE",i,compressed[i:i+6000],flush=True)
    if not any("boundaries" in m for r in result for m in r["models"].values()):
        raise RuntimeError("No trajectories propagated")

if __name__=="__main__": main()
