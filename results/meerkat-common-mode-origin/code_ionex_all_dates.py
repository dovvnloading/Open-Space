#!/usr/bin/env python3
"""Fixed 124-date independent CODE IONEX + IGRF test with 2020/21 -> 2022/23 holdout.

Use the 124 *predefined* UTC days and disjoint pulsar groups in PR #47's
dates.csv. Never select dates by RM amplitude. External independent CODE GNSS
TEC products plus IGRF-14 thin-shell RM are evaluated at every precise
observation's sky direction and (2-decimal-MJD) epoch. A model-difference
offset is evaluated on 2022/23 *without fitting to its observed RM*.

Caveat: this does NOT exactly reproduce ionFR's vertical electron profile,
geomagnetic implementation, time sampling, or source-specific 6–180 d
baselines. The result is a mechanism-sensitivity test, not a published-RM
reanalysis or a validated correction to individual pulsars.

Run with: python results/meerkat-common-mode-origin/code_ionex_all_dates.py --fetch --write
Requires numpy pandas ppigrf unlzw3, internet access to CODE.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from code_ionex_probe import (
    BASE, C_RM, LAT, LON, PRIOR, R_E, load_day, read_data, time_interp,
)

HERE = Path(__file__).resolve().parent
EPOCH = dt.datetime(1858,11,17,tzinfo=dt.timezone.utc)


def batch_ionex_rm(sample, maps, height):
    """Vectorize the expensive IGRF-14 call across pulsars in the SAME UTC day."""
    import ppigrf
    ra = np.array([z["ra_deg"] for z in sample],float)
    de = np.array([z["dec_deg"] for z in sample],float)
    mjd = np.array([z["mjd"] for z in sample],float)
    ha_hours = ((18.697374558 + 24.06570982441908*(mjd-51544.5)
                 + LON/15 - ra/15 + 12) % 24) - 12
    ha, dec, lat, lon = np.deg2rad(ha_hours*15),np.deg2rad(de),np.deg2rad(LAT),np.deg2rad(LON)
    east = -np.cos(dec)*np.sin(ha)
    north = np.sin(dec)*np.cos(lat) - np.cos(dec)*np.sin(lat)*np.cos(ha)
    up = np.sin(dec)*np.sin(lat) + np.cos(dec)*np.cos(lat)*np.cos(ha)
    altitude = np.rad2deg(np.arcsin(np.clip(up,-1,1)))
    e0 = np.array([-np.sin(lon),np.cos(lon),0])
    n0 = np.array([-np.sin(lat)*np.cos(lon),-np.sin(lat)*np.sin(lon),np.cos(lat)])
    u0 = np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    ray = east[:,None]*e0 + north[:,None]*n0 + up[:,None]*u0
    travel = -R_E*up + np.sqrt((R_E*up)**2+(R_E+height)**2-R_E**2)
    pp = R_E*u0[None,:]+travel[:,None]*ray
    pp_lat = np.rad2deg(np.arcsin(pp[:,2]/np.linalg.norm(pp,axis=1)))
    pp_lon = np.rad2deg(np.arctan2(pp[:,1],pp[:,0]))
    # One evaluation for all sightlines; ppigrf produces shape (1,N).
    reference = EPOCH + dt.timedelta(days=float(np.median(mjd)))
    br,bsouth,be = ppigrf.igrf_gc(
        np.full(len(sample),R_E+height),90-pp_lat,pp_lon,reference.replace(tzinfo=None))
    br,bsouth,be = np.ravel(br),np.ravel(bsouth),np.ravel(be)
    plat,plon=np.deg2rad(pp_lat),np.deg2rad(pp_lon)
    eu=np.stack([np.cos(plat)*np.cos(plon),np.cos(plat)*np.sin(plon),np.sin(plat)],axis=-1)
    en=np.stack([-np.sin(plat)*np.cos(plon),-np.sin(plat)*np.sin(plon),np.cos(plat)],axis=-1)
    ee=np.stack([-np.sin(plon),np.cos(plon),np.zeros(len(sample))],axis=-1)
    b_toward = -(br*np.sum(ray*eu,axis=1)-bsouth*np.sum(ray*en,axis=1)+be*np.sum(ray*ee,axis=1))
    vtec=np.array([time_interp(maps,EPOCH+dt.timedelta(days=float(t)),float(plat),float(plon))
                   for t,plat,plon in zip(mjd,pp_lat,pp_lon)])
    slant=1/np.sqrt(1-(R_E/(R_E+height)*np.cos(np.deg2rad(altitude)))**2)
    rm=C_RM*1e16*1e-9*vtec*slant*b_toward
    return rm,altitude


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--cache",type=Path,default=HERE/".ionex_cache")
    ap.add_argument("--fetch",action="store_true")
    ap.add_argument("--write",action="store_true")
    ap.add_argument("--max-days",type=int,default=124,
                    help="Debug only; any result with <124 days is NOT a full-audit claim")
    args=ap.parse_args()
    with (PRIOR/"dates.csv").open(newline="") as f:
        daily={int(float(row["mjd_utc"])):row for row in csv.DictReader(f)}
    assert len(daily)==124
    stars,baselines=read_data()
    # Group per UTC day; original PR #47 precision filters, names and hashes.
    observations=defaultdict(list)
    for n,visits in stars.items():
        if n not in baselines:
            continue
        for z in visits:
            day=int(z["mjd"])
            if day in daily and z["e"]<1.8:
                observations[day].append(z)
    receipts=[];date_rows=[]
    for day in sorted(daily)[:args.max_days]:
        date=(EPOCH+dt.timedelta(days=day)).date()
        sample=observations[day]
        (maps,height),receipt=load_day(date,args.cache,args.fetch)
        code,elevation=batch_ionex_rm(sample,maps,height)
        # Median within a named pulsar prevents repeat visits from dominating.
        byfold=[defaultdict(list),defaultdict(list)]
        skipped=0
        for z,new,alt in zip(sample,code,elevation):
            if not np.isfinite(new) or alt<10:
                skipped+=1
                continue
            fold=hashlib.sha256(z["psr"].encode()).digest()[0]%2
            byfold[fold][z["psr"]].append(float(new-z["iono"]))
        if not all(len(q)>=10 for q in byfold):
            raise ValueError(f"Too few CODE-covered distinct pulsars on {date}: {[len(q) for q in byfold]}")
        row=daily[day]
        delta=[float(np.median([np.median(v) for v in q.values()])) for q in byfold]
        date_rows.append({
            "utc_date":date.isoformat(),"mjd":day,"n_A":len(byfold[0]),"n_B":len(byfold[1]),
            "original_A":float(row["resid_A"]),"original_B":float(row["resid_B"]),
            "median_code_minus_jpl_A":delta[0],"median_code_minus_jpl_B":delta[1],
            "n_low_alt_ignored":skipped})
        receipt["date"]=date.isoformat()
        receipt["n_samples"]=len(sample)
        receipts.append(receipt)
    if len(date_rows)!=124:
        print(json.dumps({"coverage_only":True,"n":len(date_rows),
                          "note":"No 124-day inference on a partial date selection"}))
        return
    train=[z for z in date_rows if z["utc_date"] < "2022-01-01"]
    test=[z for z in date_rows if z["utc_date"] >= "2022-01-01"]
    assert len(train)==56 and len(test)==68
    # Center model differences using training days ONLY, without fitting to RM.
    centers=[float(np.median([z["median_code_minus_jpl_"+g] for z in train]))
             for g in ("A","B")]
    for z in date_rows:
        for g,cent in zip(("A","B"),centers):
            z["CODE_substituted_"+g] = (z["original_"+g]
                - (z["median_code_minus_jpl_"+g]-cent))
    def stats(v):
        a=np.array([z["original_A"] for z in v])
        b=np.array([z["original_B"] for z in v])
        ca=np.array([z["CODE_substituted_A"] for z in v])
        cb=np.array([z["CODE_substituted_B"] for z in v])
        return {"n":len(v),
                "JPL_B_rmse":float(np.sqrt(np.mean(b**2))),
                "CODE_B_rmse":float(np.sqrt(np.mean(cb**2))),
                "JPL_A_B_pearson":float(np.corrcoef(a,b)[0,1]),
                "CODE_A_B_pearson":float(np.corrcoef(ca,cb)[0,1]),
                "peer_differential_JPL_rmse":float(np.sqrt(np.mean((b-a)**2))),
                "peer_differential_CODE_rmse":float(np.sqrt(np.mean((cb-ca)**2))),
                "CODE_vs_JPL_offset_B_pearson":float(np.corrcoef(
                    b,np.array([z["median_code_minus_jpl_B"] for z in v]))[0,1]),
                "median_abs_delta_CODE_JPL_B":float(np.median(np.abs(
                    [z["median_code_minus_jpl_B"] for z in v])))}
    result={
        "description":"CODE final IONEX+IGRF independent correction on predeclared 124 source-disjoint UTC days",
        "baseline":"CODE-JPL daily difference centered at 2020-21 training median separately by source half; original PR #47 local source baselines retained; not a full ionFR re-reduction",
        "maps":len(receipts),"n_train":len(train),"n_heldout":len(test),
        "center_delta_train_A":centers[0],"center_delta_train_B":centers[1],
        "train":stats(train),"heldout":stats(test),
        "caution":"IGRF-14 thin shell and two-decimal-MJD precision; source groups share telescope; no calibrator/subarray metadata; not physical attribution. CODE-JPL delta can include geometry implementation differences.",
        "source":"https://download.aiub.unibe.ch/CODE/YYYY/...",
    }
    if args.write:
        (HERE/"code_ionex_full_checks.json").write_text(
            json.dumps(result,indent=2,sort_keys=True)+"\n")
        for filename,table in [("code_ionex_full_dates.csv",date_rows),
                               ("code_ionex_full_sources.csv",receipts)]:
            with (HERE/filename).open("w",newline="") as f:
                writer=csv.DictWriter(f,fieldnames=list(table[0]))
                writer.writeheader()
                writer.writerows(table)
    print(json.dumps({"maps":len(receipts),"train":result["train"],
                      "heldout":result["heldout"]}))


if __name__=="__main__":
    main()
