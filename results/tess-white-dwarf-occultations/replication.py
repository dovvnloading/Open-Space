"""Compare FFI triggers to the same star's independent SPOC photometry and nearby control stars."""
import argparse
import csv
import io
import json
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from astropy.io import fits
from short_pixels import acquire
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
        roster_script=(cache/f"tesscurl_sector_{s}_tp.sh").read_text()
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
                d=h[1].data
                t,tc,q,fl,bg,x,y=(d[c] for c in ("TIME","TIMECORR","QUALITY","SAP_FLUX","SAP_BKG","MOM_CENTR1","MOM_CENTR2"))
                valid=np.isfinite(t+tc)
                t,tc,q,fl,bg,x,y=(np.array(c[valid],float if i!=2 else int) for i,c in enumerate((t,tc,q,fl,bg,x,y)))
                measurements=[]
                for ev in events:
                    p=ev.get("period_days")
                    if entry["sector"]!=s and not p:continue
                    e=ev["epoch_btjd"]
                    if entry["sector"]==s:
                        sc=e-np.interp(e,times,corr)
                        shift=np.interp(sc,t-tc,tc)-np.interp(e,times,corr)
                    else:
                        # Long-baseline phase is approximate; no clock refinement in independent sectors.
                        shift=0.
                    measurements.append(dict(ffi_epoch=e,ffi_period=p,**summary(t,q,fl,bg,x,y,e,ev["duration_min"]/1440,p,shift)))
                entry["checks"]=measurements
        except Exception as exc:entry["error"]=str(exc)
        result.append(entry)
    (out/(tic+".json")).write_text(json.dumps(dict(tic=tic,gaia=a["source_id"],sector=s,comparison_roster=meta,products=result),allow_nan=False)+"\n")
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("curves","parent","inputs","inventory","cache","out"):p.add_argument(k,type=Path)
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True);a.cache.mkdir(exist_ok=True,parents=True)
    stars={r["GaiaEDR3"]:r for r in csv.DictReader(a.parent.open())};roster=list(csv.DictReader(a.inventory.open()))
    jobs=[]
    for f in sorted(a.curves.glob("*.json")):
        row=json.loads(f.read_text())
        if row.get("status")!="included" or not(row["search"]["periodic"] or row["search"]["isolated"]):continue
        if (a.out/(row["tic"]+".json")).exists():continue
        star=stars[row["source_id"]];row.update(ra_obs=float(star["RA_ICRS"]),dec_obs=float(star["DE_ICRS"]));jobs.append(row)
    with ThreadPoolExecutor(4) as ex:list(ex.map(lambda r:check_host(r,a.inputs,a.cache,roster,a.out),jobs))
