"""Refit one asteroid-like slow object and two marginal CODE comet entries.

This replay keeps the original nine-object input set intact. Freeze each newly
acquired measurement set and report failures without replacing them with catalogue orbits.
"""
import base64
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import urllib.parse

import numpy as np
from jplephem.spk import SPK

ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
import acquire
import refit as fit

OUT = ROOT / "evidence"
OUT.mkdir(exist_ok=True)
CACHE = ROOT / "cache"
CACHE.mkdir(exist_ok=True)
acquire.CACHE=fit.CACHE=CACHE
fit.DATA=OUT
TARGETS = {"2024-U2":("A/2024 U2","K24U020","54512605"),
           "2011-M1":("C/2011 M1","K11M010","1003135"),
           "2020-W5":("C/2020 W5","K20W050","1003754")}
BULK="https://www.minorplanetcenter.net/iau/ECS/MPCAT-OBS/CmtObs.txt.gz"
MPC="https://data.minorplanetcenter.net/api/get-obs"
PLANETS="https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440.bsp"
STATIONS="https://minorplanetcenter.net/iau/lists/ObsCodes.html"
KERNEL_SHA="a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"

def archive(name,b):
    encoded=base64.b64encode(gzip.compress(b,mtime=0)).decode()
    for start in range(0,len(encoded),6000):
        print("SLOW_CHUNK",name,start,encoded[start:start+6000],flush=True)

def observations(obj,name,packed,receipts):
    body=json.dumps({"desigs":[name],"output_format":["OBS80"]}).encode()
    source="live MPC Observations API"
    try:
        wire=acquire.retrieve(MPC,obj+".json",receipts,body)
        field=json.loads(wire)[0]["OBS80"]
        lines=acquire.select(field.splitlines() if isinstance(field,str) else field,[packed])
    except (ValueError,OSError,TypeError,KeyError,IndexError):
        lines=[]
    if not lines:
        wire=acquire.retrieve(BULK,"bulk.gz",receipts)
        lines=acquire.select(gzip.decompress(wire).decode("ascii","replace").splitlines(),[packed])
        source="MPC CmtObs bulk (fallback)"
    if not lines:
        raise ValueError("no MPC astrometry with matching packed designation")
    b=("\n".join(lines)+"\n").encode()
    (OUT/(obj+".obs")).write_bytes(b)
    archive(obj+".obs.gz",b)
    return {"source":source,"sha256":hashlib.sha256(b).hexdigest(),"records":len(lines),
            "start_date":min(l[15:32] for l in lines).strip(),"end_date":max(l[15:32] for l in lines).strip()}

def initial_state(obj,spk,epoch,receipts):
    params={"format":"json","COMMAND":"'"+spk+"'","OBJ_DATA":"NO","EPHEM_TYPE":"VECTORS",
            "CENTER":"'500@0'","TLIST":f"'{epoch:.8f}'","CSV_FORMAT":"YES","REF_PLANE":"FRAME",
            "VEC_TABLE":"2","OUT_UNITS":"AU-D","VEC_CORR":"NONE"}
    url="https://ssd.jpl.nasa.gov/api/horizons.api?"+urllib.parse.urlencode(params)
    wire=acquire.retrieve(url,obj+"-seed.json",receipts)
    text=json.loads(wire)["result"]
    row=next(csv.reader(io.StringIO(text.split("$$SOE")[1].split("$$EOE")[0].strip())))
    return np.array([float(x) for x in row[2:8]])

def measure(obj,name,packed,spk,eph,receipts,entry):
    entry["input"]=observations(obj,name,packed,receipts)
    t,observers,east,north,sig,groups,rejected=fit.observations(obj,eph)
    epoch=float(np.ceil(max(t))+1)
    seed=initial_state(obj,spk,epoch,receipts)
    positions=fit.ephemeris(eph,min(t)-fit.OBS_TIME_LIMIT-2,epoch+1,.25)
    entry.update(supported=len(t),n_station_nights=len(set(groups)),selection_rejections=rejected,
                 fit_epoch_jd_tdb=epoch,initial_state=seed.tolist(),models={})
    models=["gravity"] if name.startswith("A/") else ["gravity","radial_r2","marsden_outward"]
    for kind in models:
        record={}
        entry["models"][kind]=record
        try:
            nng={"gravity":0,"radial_r2":1,"marsden_outward":2}[kind]
            x,cov,mask,res,redchi,pred=fit.fit(np.r_[seed,np.zeros(nng)],epoch,t,observers,east,north,
                    sig,positions,outward=kind=="marsden_outward",return_predictor=True)
            record.update(state_ssb_icrf=x.tolist(),covariance=cov.tolist(),used=int(mask.sum()),
                          rms_as=float(np.sqrt(np.mean(res[mask]**2))),reduced_chi2=float(redchi))
            if kind=="gravity": seed=x.copy()
            record["boundaries"]=fit.inbound(obj,x,cov,epoch,eph)
            boundary=record["boundaries"].get("250au",{})
            if boundary.get("inverse_a_au_inv",0)<0 and kind=="gravity":
                record["holdouts"]=fit.check_holdouts(obj,x,epoch,eph,t,observers,east,north,sig,groups,positions)
            if boundary.get("inverse_a_au_inv",0)<0:
                try:
                    record["zero_energy_profile"]=fit.parabolic_profile(epoch,eph,x,kind,mask,sig,pred)
                except Exception as exc: record["zero_energy_profile"]={"failed":str(exc)}
        except Exception as exc:
            record["failed"]=str(exc)

def main():
    receipts=[]; rows={}
    kernel=acquire.retrieve(PLANETS,"de440.bsp",receipts)
    if hashlib.sha256(kernel).hexdigest()!=KERNEL_SHA: raise ValueError("DE440 checksum mismatch")
    sites=acquire.retrieve(STATIONS,"stations.txt",receipts)
    if b"Greenwich" not in sites: raise ValueError("missing observatory table")
    with SPK.open(str(CACHE/"de440.bsp")) as eph:
        for obj,(name,packed,spk) in TARGETS.items():
            rows[obj]={"name":name,"spkid":spk}
            try: measure(obj,name,packed,spk,eph,receipts,rows[obj])
            except Exception as exc: rows[obj]["failed"]=str(exc)
            print("SLOW_RESULT",obj,json.dumps(rows[obj],allow_nan=False),flush=True)
    output={"objects":rows,"queries":receipts}
    b=(json.dumps(output,indent=2,allow_nan=False)+"\n").encode()
    (OUT/"results.json").write_bytes(b)
    archive("results.json.gz",b)
    print("SLOW_QUERIES",json.dumps(receipts),flush=True)

if __name__=="__main__": main()
