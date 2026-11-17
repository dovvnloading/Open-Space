"""Compare FFI triggers to the same star's independent SPOC photometry and nearby control stars."""
import argparse
import csv
import hashlib
import gzip
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from astropy.io import fits
from short_pixels import acquire
from ffi_acquire import download, propagated
from ffi_images import event_masks

def direct(name):
    sector,sid=name.split("-")[1:3]
    return f"https://stpubdata.s3.amazonaws.com/tess/public/tid/{sector}/"+"/".join(sid[i:i+4] for i in (0,4,8,12))+"/"+name
def summary(t,qual,flux,sky,cx,cy,center,dur,p,shift):
    mid,l,r=event_masks(t,center+shift,dur,p)
    q=(qual==0)&np.isfinite(flux)
    raw=mid.sum()
    inside=q&mid;left=q&l;right=q&r
    valid=int(inside.sum())
    out=dict(raw=int(raw),quality_zero_finite=valid,n_left=int(left.sum()),n_right=int(right.sum()),center_shift_sec=round(shift*86400,3))
    if raw:
        a,b=np.unique(qual[mid],return_counts=True);out["quality_counts"]={str(int(k)):int(n) for k,n in zip(a,b)}
    if not (valid and left.any() and right.any()):return out
    # Balanced shoulders on either side to reject a one-sided step; record sides too.
    inn,bef,aft=(float(np.mean(flux[m])) for m in (inside,left,right))
    baseline=(bef+aft)/2
    def change(x):
        if any(not np.isfinite(x[m]).any() for m in (inside,left,right)):return None
        return float(np.nanmean(x[inside])-.5*(np.nanmean(x[left])+np.nanmean(x[right])))
    if baseline>0:
        out.update(relative_depth=1-inn/baseline,left_right_fraction=(bef-aft)/baseline,sky_shift_e_s=change(sky),centroid_col_shift=change(cx),centroid_row_shift=change(cy))
    if not p:
        context=q&(abs(t-center-shift)>2*dur)&(abs(t-center-shift)<max(.1,4*dur))
        if context.sum()>5:
            ref=np.median(flux[context])
            if ref>0:out.update(center_vs_local_med=float(inn/ref-1),shoulders_vs_local_med=float(baseline/ref-1))
    return out
def check_host(a,inputs,cache,inventory,out):
    events=[dict(kind=k,**e) for k in ("periodic","isolated") for e in a["search"][k]]
    if not events:return
    tic,s=a["tic"],a["sector"]
    ffi=inputs/f"{a['source_id']}-s{s}.fits"
    with fits.open(ffi) as h:
        times=h[1].data["TIME"].copy();corr=h[1].data["TIMECORR"].copy()
        good=np.isfinite(times+corr);times,corr=times[good],corr[good]
    own=[r for r in inventory if r["tic"]==tic and (int(r["sector"])==s or any(e["kind"]=="periodic" for e in events))]
    # Add nearest three two-minute targets on the same camera and CCD, excluding the WD.
    roster_url=f"https://tess.mit.edu/public/target_lists/2m/all_targets_S{s:03d}_v1.csv"
    roster=cache/f"sector{s}-mit.csv";meta={"url":roster_url}
    try:
        meta["sha256"]=acquire(roster,roster_url)
        text=[r.strip() for r in roster.read_text().splitlines() if r.strip() and not r.startswith("#")]
        targets=[]
        for line in text:
            f=[x.strip() for x in line.split(",")]
            if not f[0].isdigit() or len(f)<6:continue
            if int(f[1])==a["camera"] and int(f[2])==a["ccd"] and f[0]!=tic:
                dist=3600*np.hypot((float(f[4])-a["ra_obs"])*np.cos(np.radians(a["dec_obs"])),float(f[5])-a["dec_obs"])
                if dist>120:targets.append((dist,f[0]))
        script=cache/f"tesscurl_sector_{s}_tp.sh"
        meta["mast_roster_sha256"]=hashlib.sha256(script.read_bytes()).hexdigest()
        roster_script=script.read_text()
        for dist,other in sorted(targets)[:3]:
            m=re.search(r"-o (tess\S+-"+other.zfill(16)+r"-\S+_tp\.fits) ",roster_script)
            if m:own.append(dict(tic=other,sector=s,cadence_s=120,filename=m[1],separation_arcsec=float(dist)))
        meta["comparison_products"]=len([r for r in own if r["tic"]!=tic])
    except Exception as exc:meta["error"]=str(exc)
    result=[]
    for row in own:
        name=row["filename"].replace("_fast-tp.fits","_fast-lc.fits").replace("_tp.fits","_lc.fits")
        url=direct(name);path=cache/name
        entry=dict(filename=name,url=url,tic=row["tic"],sector=int(row["sector"]),cadence_s=int(row["cadence_s"]),separation_arcsec=row.get("separation_arcsec",0))
        try:
            entry["sha256"]=acquire(path,url)
            with fits.open(path) as h:
                if (int(h[0].header["TICID"])!=int(row["tic"]) or
                        int(h[0].header["SECTOR"])!=int(row["sector"]) or
                        abs(86400*h[1].header["TIMEDEL"]-int(row["cadence_s"]))>0.1):
                    raise ValueError("SPOC product identity or cadence mismatch")
                entry["spoc_aperture_target_fraction"]=h[1].header.get("CROWDSAP")
                d=h[1].data
                t,tc,q,fl,bg,x,y=(d[c] for c in ("TIME","TIMECORR","QUALITY","SAP_FLUX","SAP_BKG","MOM_CENTR1","MOM_CENTR2"))
                valid=np.isfinite(t+tc)
                t,tc,q,fl,bg,x,y=(np.array(c[valid],float if i!=2 else int) for i,c in enumerate((t,tc,q,fl,bg,x,y)))
                # Use the FFI barycentric time system for every simultaneous exposure.
                # The interpolated correction varies over a visit and between target lines of sight.
                if entry["sector"]==s:
                    spacecraft=t-tc
                    mapped_time=spacecraft+np.interp(spacecraft,times-corr,corr)
                else:
                    mapped_time=t
                measurements=[]
                for ev in events:
                    p=ev.get("period_days")
                    if entry["sector"]!=s and not p:continue
                    e=ev["epoch_btjd"]
                    dur=ev["duration_min"]/1440
                    check=summary(mapped_time,q,fl,bg,x,y,e,dur,p,0.)
                    if p and entry["sector"]==s and row["tic"]==tic:
                        # Timestamp-specific quality and the two local flanks for each predicted orbit.
                        k=np.arange(np.ceil((mapped_time.min()-e)/p),np.floor((mapped_time.max()-e)/p)+1)
                        cycles=[]
                        for n in k:
                            c=e+n*p
                            sl=slice(np.searchsorted(mapped_time,c-4.1*max(dur,.03)),
                                     np.searchsorted(mapped_time,c+4.1*max(dur,.03)))
                            cycles.append(dict(cycle=int(n),**summary(mapped_time[sl],q[sl],fl[sl],
                                bg[sl],x[sl],y[sl],c,dur,None,0.)))
                        check["cycles"]=cycles
                    measurements.append(dict(ffi_epoch=e,ffi_period=p,**check))
                entry["checks"]=measurements
        except Exception as exc:entry["error"]=str(exc)
        result.append(entry)
    (out/(tic+".json")).write_text(json.dumps(dict(tic=tic,gaia=a["source_id"],sector=s,comparison_roster=meta,products=result),allow_nan=False)+"\n")
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("curves","parent","inputs","inventory","cache","out"):p.add_argument(k,type=Path)
    p.add_argument("--receipts",type=Path,help="Archive receipt list to materialize the original FFI stamps")
    p.add_argument("--bulk",type=Path,help="Archive roster receipt list to materialize MAST neighbour rosters")
    p.add_argument("--pack",type=Path,help="Package completed records for this input cohort in TIC order")
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True);a.cache.mkdir(exist_ok=True,parents=True)
    a.inputs.mkdir(exist_ok=True,parents=True)
    stars={r["GaiaEDR3"]:r for r in csv.DictReader(a.parent.open())};roster=list(csv.DictReader(a.inventory.open()))
    rows=json.loads(a.curves.read_text()) if a.curves.is_file() else [json.loads(f.read_text()) for f in sorted(a.curves.glob("*.json"))]
    jobs=[]
    for row in rows:
        if row.get("status")!="included" or not(row["search"]["periodic"] or row["search"]["isolated"]):continue
        if (a.out/(row["tic"]+".json")).exists():continue
        star=stars[row["source_id"]];ra,dec=propagated(star,row["coordinate_epoch"])
        row.update(ra_obs=ra,dec_obs=dec);jobs.append(row)
    if a.bulk:
        for meta in json.loads(a.bulk.read_text())["scripts"]:
            if meta["sector"] in {r["sector"] for r in jobs} and meta["cadence_s"]==120:
                script=a.cache/Path(meta["url"]).name
                if acquire(script,meta["url"])!=meta["sha256"]:raise ValueError("MAST roster changed: "+script.name)
    if a.receipts:
        refs={r["file"]:r for r in json.loads(a.receipts.read_text())}
        def get_stamp(r):
            path=a.inputs/f"{r['source_id']}-s{r['sector']}.fits"
            ref=refs[path.name+".json"]
            if not path.exists():download(path,ref["url"])
            with fits.open(path) as h:
                for key,digest in r["arrays_sha256"].items():
                    if hashlib.sha256(h[1].data[key].tobytes()).hexdigest()!=digest:
                        raise ValueError("Frozen FITS science array changed: "+path.name+" "+key)
        with ThreadPoolExecutor(4) as ex:list(ex.map(get_stamp,jobs))
    with ThreadPoolExecutor(4) as ex:list(ex.map(lambda r:check_host(r,a.inputs,a.cache,roster,a.out),jobs))
    if a.pack:
        hosts=sorted({r["tic"] for r in rows if r.get("status")=="included"
                      and (r["search"]["periodic"] or r["search"]["isolated"])})
        records=[json.loads((a.out/(tic+".json")).read_text()) for tic in hosts]
        a.pack.write_bytes(gzip.compress(json.dumps(records,separators=(",",":"),allow_nan=False).encode(),mtime=0))
