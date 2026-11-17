"""Gravity-screen the disjoint well-observed portion of the frozen e >= 1 catalogue."""
import argparse
import base64
from collections import defaultdict
import concurrent.futures
import gzip
import hashlib
import importlib.util
import json
import multiprocessing
import os
from pathlib import Path
import signal
import sys
import time

import numpy as np
from jplephem.spk import SPK

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/"slow"))
from replay import acquire, conic_seed, fit, KERNEL_SHA, PLANETS, STATIONS
BATCH=ROOT.parents[1]/"batch"
spec=importlib.util.spec_from_file_location("baseline_batch",BATCH/"refit.py")
batch=importlib.util.module_from_spec(spec)
spec.loader.exec_module(batch)
CACHE=ROOT/"cache"
OUT=ROOT/"evidence"
CACHE.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)
acquire.CACHE=fit.CACHE=CACHE
fit.DATA=CACHE
SUN_MU_KM=132712440041.93938/149597870.7
SPECIFIC={"1877 G2","1997 P2","2008 J4","1980 E1","2024 L5","2025 V1","2011 M1","2020 W5"}
ORIGINAL_EPHEMERIS=fit.ephemeris
PAST=None
MAX_POSITIONS=96

def digest(blob):
    return hashlib.sha256(blob).hexdigest()

def limited_ephemeris(eph,first,last,step):
    if step==1 and PAST is not None and first>=PAST.x[0] and last<=PAST.x[-1]:
        return PAST
    return ORIGINAL_EPHEMERIS(eph,first,last,step)

def representatives(lines,sites,short=False):
    groups=defaultdict(list)
    for raw in lines:
        line=raw.decode("ascii","replace").rstrip("\r\n")
        if len(line)>=80 and line[14] in " CBAOPTMc" and line[77:80] in sites:
            # Reject the continuation of two-line records and invalid dates.
            if line[15:19].isdigit() and line[15:25].count(" ")==2:
                groups[(line[15:25],line[77:80])].append(line[:80])
    reps=[]
    for _,group in sorted(groups.items()):
        group=sorted(group,key=lambda s:(s[15:32],s[32:56]) if short else (s[15:32],))
        if short:
            # Preserve angular-rate information within sparse nights.
            seen=set()
            raw=[]
            for line in group:
                if line[15:56] not in seen:
                    seen.add(line[15:56])
                    raw.append(line)
            k=min(8,len(raw))
            reps.extend([raw[round(i*(len(raw)-1)/(k-1))] for i in range(k)] if k>1 else raw)
        else:
            reps.append(group[(len(group)-1)//2])
    if len(reps)>MAX_POSITIONS:
        reps=[reps[round(i*(len(reps)-1)/(MAX_POSITIONS-1))] for i in range(MAX_POSITIONS)]
    return len(groups),reps

def abort(signum,frame):
    raise TimeoutError("180 s per-object wall-clock budget")

def measure(item):
    row,lines,short=item
    key=str(row["spkid"])
    names=["q","e","i","om","w","tp","epoch"]
    start=time.monotonic()
    rec={"spkid":key,"name":row["full_name"].strip(),"reported_arc_days":float(row["data_arc"]) if row["data_arc"] else None,
         "bulk_lines":len(lines),"bulk_astrometry_sha256":digest(b"".join(lines)) if lines else None}
    signal.signal(signal.SIGALRM,abort)
    signal.alarm(180)
    minimum=6 if short else 10
    span_min=2 if short else 30
    try:
        sites=fit.stations()
        nights,selected=representatives(lines,sites,short)
        raw=("\n".join(selected)+"\n").encode("ascii")
        rec.update(available_station_nights=nights,selected_mpc80=selected,
                   selected_sha256=digest(raw) if selected else None)
        if not short and nights<10: raise ValueError("fewer than ten ground optical station-nights")
        if not lines: raise ValueError("no matching packed-designation records in CmtObs snapshot")
        if not selected: raise ValueError("no supported ground optical station-nights")
        if nights<(3 if short else 10): raise ValueError("insufficient ground optical station-nights")
        if len(selected)<minimum: raise ValueError("insufficient distinct ground optical positions")
        if not all(row.get(k) is not None for k in names):
            raise ValueError("frozen SBDB conic missing fields")
        (CACHE/(key+".obs")).write_bytes(raw)
        with SPK.open(str(CACHE/"de440.bsp")) as eph:
            t,observer,east,north,errors,groups,dropped=fit.observations(key,eph,minimum_positions=6 if short else 7)
            epoch=float(np.ceil(t.max())+1)
            span=float(np.ptp(t))
            rec.update(supported=len(t),supported_arc_days=span,parser_dropped=dropped)
            if not short and (len(t)<10 or span<30):
                raise ValueError("supported arc shorter than 30 days or ten positions")
            if len(t)<minimum or span<span_min or (short and len(set(groups))<3):
                raise ValueError("supported arc below required span, positions, or distinct station-nights")
            el={k:float(row[k]) for k in names}
            initial=conic_seed(el,el["epoch"],epoch,eph)
            eph_short=fit.ephemeris(eph,t.min()-fit.OBS_TIME_LIMIT-2,epoch+1,.25)
            state,cov,keep,res,chi=fit.fit(initial,epoch,t,observer,east,north,errors,eph_short,minimum_positions=6 if short else None)
            # Check for the antipodal tangent-plane degeneracy.
            sol=fit.integrate(state,epoch,t.min()-fit.OBS_TIME_LIMIT,eph_short).sol
            emit=t.copy()
            for _ in range(4):
                d=sol(emit)[:3].T-observer
                emit=t-np.linalg.norm(d,axis=1)/fit.C
            d=sol(emit)[:3].T-observer
            alignment=np.einsum("ij,ij->i",d/np.linalg.norm(d,axis=1)[:,None],np.cross(east,north))
            if np.min(alignment[keep])<.99: raise ValueError("direction check failed")
            rec.update(used=int(keep.sum()),used_arc_days=float(np.ptp(t[keep])),
                       used_count_sites=len({g[0] for g,flag in zip(groups,keep) if flag}),
                       excluded_indices=np.where(~keep)[0].tolist(),rms_coordinate_arcsec=float(np.sqrt(np.mean(res[keep]**2))),
                       reduced_chi2=float(chi),fit_epoch_jd_tdb=epoch,state_ssb_icrf=state.tolist(),
                       covariance=cov.tolist())
            used_groups={g for g,flag in zip(groups,keep) if flag}
            if not short and (keep.sum()<10 or np.ptp(t[keep])<30 or np.ptp(t[keep])<.7*span):
                raise ValueError("retained fit spans <30 d, <70% of arc, or <10 points")
            if keep.sum()<minimum or (short and len(used_groups)<3) or np.ptp(t[keep])<span_min or np.ptp(t[keep])<.7*span:
                raise ValueError("retained arc below span, 70% coverage, positions, or station-night floor")
            rec["boundaries"]=fit.inbound(key,state,cov,epoch,eph)
            b=rec["boundaries"].get("250au")
            rec["status"]="inbound 250 au" if b and b["inbound_radial_au_d"]<0 else "no inbound 250 au"
            if short and b and (not np.isfinite(b["formal_sigma_inverse_a"]) or b["formal_sigma_inverse_a"]<=0):
                raise ValueError("nonpositive or nonfinite energy uncertainty")
    except Exception as exc:
        rec.update(status="unresolved",reason=str(exc))
    finally:
        signal.alarm(0)
    rec["elapsed_seconds"]=round(time.monotonic()-start,3)
    return rec

def main(src,shard,parts,short=False):
    global PAST
    rows,url,cat_sha=batch.cat()
    source=src/"CmtObs.txt.gz"
    receipt=json.loads((src/"source.json").read_text())
    if digest(source.read_bytes())!=receipt["sha256_compressed"]:
        raise ValueError("CmtObs snapshot hash mismatch")
    live=BATCH/"evidence/latest_mpc_refits.json.gz"
    recent={str(r["spkid"]) for r in json.loads(gzip.decompress(live.read_bytes()))}
    keys=defaultdict(set)
    for r in rows:
        for alias in batch.aliases(r): keys[alias].add(str(r["spkid"]))
    eligible=[]
    for r in rows:
        speed=np.sqrt(SUN_MU_KM*(float(r["e"])-1)/float(r["q"]))
        long_arc=float(r["data_arc"] or -1)>=30
        if speed<1 and long_arc!=short and str(r["spkid"]) not in recent and r["pdes"] not in SPECIFIC:
            eligible.append(r)
    assigned=sorted(eligible,key=lambda r:int(r["spkid"]))[shard::parts]
    ids={str(r["spkid"]) for r in assigned}
    groups=defaultdict(list)
    with gzip.open(source,"rb") as f:
        for line in f:
            aliases=keys.get(line[5:12].decode("ascii","replace").strip() or line[:5].decode("ascii","replace").strip(),set())
            if len(aliases)==1:
                ident=next(iter(aliases))
                if ident in ids: groups[ident].append(line)
    queries=[]
    if digest(acquire.retrieve(PLANETS,"de440.bsp",queries))!=KERNEL_SHA:
        raise ValueError("DE440 hash mismatch")
    if b"Greenwich" not in acquire.retrieve(STATIONS,"stations.txt",queries):
        raise ValueError("invalid observatory list")
    # Construct the common one-day DE440 spline once before forking. It has the
    # same integer-day knots as refit.inbound. It is read-only in all workers.
    with SPK.open(str(CACHE/"de440.bsp")) as eph:
        lo=max(seg.start_jd for seg in eph.segments if seg.center==0)+5
        hi=max(float(r["epoch"] or 0) for r in eligible)+50
        hi=max(hi,2461500.)
        hi=min(hi,min(seg.end_jd for seg in eph.segments if seg.center==0)-4)
        PAST=ORIGINAL_EPHEMERIS(eph,lo,hi,1.)
    fit.ephemeris=limited_ephemeris
    receipt.update(catalogue_sha256=cat_sha,catalogue_url=url,latest_exclusion_archive_sha256=digest(live.read_bytes()),
                   total_eligible=len(eligible),shard=shard,parts=parts,assigned=len(assigned),additional_queries=queries)
    print("CENSUS_SOURCE",json.dumps(receipt),flush=True)
    dest=OUT/f"shard-{shard:02}.json.gz"
    records=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context("fork")) as pool:
        futures=[pool.submit(measure,(r,groups[str(r["spkid"])],short)) for r in assigned]
        for future in concurrent.futures.as_completed(futures):
            r=future.result()
            records.append(r)
            boundary=r.get("boundaries",{}).get("250au",{})
            print("CENSUS_ROW",r["spkid"],r["status"],r.get("reason",""),boundary.get("inverse_a_au_inv"),boundary.get("formal_sigma_inverse_a"),r["elapsed_seconds"],flush=True)
            records.sort(key=lambda v:int(v["spkid"]))
            dest.write_bytes(gzip.compress(json.dumps({"receipt":receipt,"objects":records},separators=(",",":"),allow_nan=False).encode(),mtime=0))
    raw=base64.b64encode(dest.read_bytes()).decode()
    for i in range(0,len(raw),5000):
        print("CENSUS_DATA",i,raw[i:i+5000],flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("snapshot",type=Path)
    p.add_argument("shard",type=int)
    p.add_argument("parts",type=int)
    args=p.parse_args()
    main(args.snapshot,args.shard,args.parts)
