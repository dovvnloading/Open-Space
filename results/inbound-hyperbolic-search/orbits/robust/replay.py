"""Refit formal negative-energy flags using a denser MPC slice and activity models."""
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import base64
import gzip
import hashlib
import importlib.util
import json
import multiprocessing
from pathlib import Path
import signal
import sys

import numpy as np
from jplephem.spk import SPK

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("ground_census",ROOT.parent/"census/replay.py")
ground=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=ground
spec.loader.exec_module(ground)
acquire,fit=ground.acquire,ground.fit
CACHE=ROOT/"cache"
OUT=ROOT/"evidence"
CACHE.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)
fit.CACHE=acquire.CACHE=fit.DATA=CACHE

def digest(b): return hashlib.sha256(b).hexdigest()

def dense_selection(rows,sites):
    group=defaultdict(dict)
    for l in rows:
        if l[14] not in " CBAOPTMc" or l[77:80] not in sites or not l[15:19].isdigit(): continue
        group[(l[15:25],l[77:80])][l[15:56]]=l
    all_positions=sorted((l for g in group.values() for l in g.values()),key=lambda l:(l[15:32],l[77:80],l[32:56]))
    if len(all_positions)<=1200: return len(all_positions),all_positions
    rows=[]
    for k,g in sorted(group.items()):
        ordered=sorted(g.values(),key=lambda l:l[15:32])
        n=min(len(ordered),8)
        rows.extend(ordered if n<2 else [ordered[round(i*(len(ordered)-1)/(n-1))] for i in range(n)])
    if len(rows)>1200: rows=[rows[round(i*(len(rows)-1)/1199)] for i in range(1200)]
    return len(all_positions),rows

def deadline(_signum,_frame): raise TimeoutError("300 s per-object model budget")

def one(item):
    screen,row,raw=item
    spk=str(row["spkid"])
    key="check-"+spk
    rec={"spkid":spk,"name":row["full_name"].strip(),"screen":screen,"models":{}}
    signal.signal(signal.SIGALRM,deadline)
    signal.alarm(300)
    try:
        available,obs=dense_selection(raw,fit.stations())
        b=("\n".join(obs)+"\n").encode("ascii")
        (CACHE/(key+".obs")).write_bytes(b)
        rec.update(available_ground_positions=available,selected_mpc80=obs,selected_sha256=digest(b),selected=len(obs))
        with SPK.open(str(CACHE/"de440.bsp")) as eph:
            t,observer,east,north,sigma,groups,dropped=fit.observations(key,eph,minimum_positions=6)
            epoch=float(np.ceil(t.max())+1)
            el={k:float(row[k]) for k in ["q","e","i","om","w","tp","epoch"]}
            seed=ground.conic_seed(el,el["epoch"],epoch,eph)
            positions=fit.ephemeris(eph,t.min()-fit.OBS_TIME_LIMIT-2,epoch+1,.25)
            rec.update(supported=len(t),arc_days=float(np.ptp(t)),preperihelion_count=int(np.sum(t<el["tp"])),
                       fit_epoch_jd_tdb=epoch,stations=len({g[0] for g in groups}),station_nights=len(set(groups)))
            jobs=[("gravity",0,None),("radial_r2",1,None),("marsden_outward",2,None)]
            if (t<el["tp"]).sum()>=12 and np.ptp(t[t<el["tp"]])>=30 and (t>=el["tp"]).sum()>=3:
                jobs.append(("preperihelion_gravity",0,t<el["tp"]))
            if t.max()<2436934.5:
                jobs.append(("historical_10arcsec",0,None))
            for kind,ng,under in jobs:
                part=np.ones(len(t),bool) if under is None else under
                rs={}
                rec["models"][kind]=rs
                try:
                    weight=sigma[part]*(2 if kind=="historical_10arcsec" else 1)
                    x,cov,keep,res,chi,pred=fit.fit(np.r_[seed,np.zeros(ng)],epoch,t[part],observer[part],east[part],north[part],
                                            weight,positions,outward=ng==2,return_predictor=True,
                                            minimum_positions=6 if ng==0 else None)
                    if kind=="gravity": seed=x.copy()
                    rs.update(used=int(keep.sum()),rms_coordinate_arcsec=float(np.sqrt(np.mean(res[keep]**2))),
                              reduced_chi2=float(chi),state_ssb_icrf=x.tolist(),covariance=cov.tolist(),
                              used_arc_days=float(np.ptp(t[part][keep])),excluded_indices_wrt_supported=np.where(part)[0][~keep].tolist())
                    rs["boundaries"]=fit.inbound(key,x,cov,epoch,eph)
                    q=rs["boundaries"].get("250au") or {}
                    if kind=="gravity":
                        rs["holdouts"]=fit.check_holdouts(key,x,epoch,eph,t,observer,east,north,sigma,groups,positions)
                    if q and q["inverse_a_au_inv"]< -2.5*q["formal_sigma_inverse_a"] and kind=="marsden_outward":
                        try: rs["parabolic_profile"]=fit.parabolic_profile(epoch,eph,x,kind,keep,sigma,pred)
                        except Exception as exc: rs["profile_failed"]=str(exc)
                except TimeoutError: raise
                except Exception as exc: rs["failed"]=str(exc)
    except Exception as exc:
        rec["error"]=str(exc)
    finally: signal.alarm(0)
    return rec

def main(src,index,parts):
    receipt=json.loads((src/"source.json").read_text())
    if digest((src/"CmtObs.txt.gz").read_bytes())!=receipt["sha256_compressed"]: raise ValueError("bulk hash mismatch")
    rows,_,catsha=ground.batch.cat()
    by_id={str(r["spkid"]):r for r in rows}
    aliases=defaultdict(set)
    for r in rows:
        for k in ground.batch.aliases(r): aliases[k].add(str(r["spkid"]))
    queue=(ROOT/"current.txt").read_text().strip()
    if queue not in {"first","second"}: raise ValueError("invalid target list name")
    wire=(ROOT/("targets-"+queue+".json")).read_bytes()
    target=json.loads(wire)[index::parts]
    ids={str(r["spkid"]) for r in target}
    recs=defaultdict(list)
    with gzip.open(src/"CmtObs.txt.gz","rt",encoding="ascii",errors="replace") as inp:
        for line in inp:
            who=aliases.get(line[5:12].strip() or line[:5].strip(),set())
            if len(who)==1:
                spk=next(iter(who))
                if spk in ids and len(line)>=80: recs[spk].append(line[:80])
    downloads=[]
    if digest(acquire.retrieve(ground.PLANETS,"de440.bsp",downloads))!=ground.KERNEL_SHA: raise ValueError("DE440 hash mismatch")
    if b"Greenwich" not in acquire.retrieve(ground.STATIONS,"stations.txt",downloads): raise ValueError("bad station file")
    inputs=[(r,by_id[r["spkid"]],recs[r["spkid"]]) for r in target]
    out={"catalogue_sha256":catsha,"source":receipt,"queries":downloads,"target_list_sha256":digest(wire),"objects":[]}
    dest=OUT/f"{queue}-part-{index}.json.gz"
    with ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context("fork")) as pool:
        for row in pool.map(one,inputs):
            out["objects"].append(row)
            print("CHECKED",row["spkid"],[(k,m.get("boundaries",{}).get("250au",{}).get("inverse_a_au_inv"),m.get("failed")) for k,m in row["models"].items()],flush=True)
            dest.write_bytes(gzip.compress(json.dumps(out,separators=(",",":"),allow_nan=False).encode(),mtime=0))
    if not dest.exists(): dest.write_bytes(gzip.compress(json.dumps(out).encode(),mtime=0))
    b=base64.b64encode(dest.read_bytes()).decode()
    for i in range(0,len(b),5000): print("ROBUST_DATA",i,b[i:i+5000],flush=True)

if __name__=="__main__": main(Path(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]))
