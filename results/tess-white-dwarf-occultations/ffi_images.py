"""Remeasure threshold events and the strongest individual subthreshold dip using the pixels."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

def module(path):
    spec=importlib.util.spec_from_file_location("wd_search",path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def event_masks(t,epoch,duration,period):
    d=(t-epoch+.5*period)%period-.5*period if period else (t-epoch)
    inside=np.abs(d)<=duration/2+1e-8
    left=(d>=-1.8*duration)&(d<=-.8*duration)
    right=(d>=.8*duration)&(d<=1.8*duration)
    return inside,left,right

def describe(rec,search_code,inputs,curves):
    tic=rec["tic"]
    source=rec["source_id"]
    sector=rec["sector"]
    with np.load(curves/rec["lc_file"]) as lc:
        t,f=lc["time"],lc["flux"]
        row,ap,cx,cy=lc["pixel_row"],lc["aperture"],lc["centroid_x"],lc["centroid_y"]
        dt=float(lc["cadence_s"])/86400
    report=search_code.search(t,f,dt*86400)
    base=search_code.flatten(t,f)
    ok=np.isfinite(base)
    t,f,row,cx,cy,base=[x[ok] for x in (t,f,row,cx,cy,base)]
    resid=f/base-1.
    strongest=None
    for width in np.unique(np.maximum(dt,np.asarray(search_code.DURATIONS_MIN)/1440)):
        depth,z=search_code.box_scores(t,resid,dt,width,search_code.mad(resid))
        if not np.isfinite(z).any(): continue
        i=np.nanargmax(z)
        if strongest is None or z[i]>strongest["score"]:
            strongest=dict(epoch_btjd=float(t[i]),duration_min=width*1440,
                           score=float(z[i]),depth=float(depth[i]),diagnostic="strongest single box")
    events=[dict(ev,diagnostic="periodic trigger") for ev in report["periodic"]]
    events.extend(dict(ev,diagnostic="individual trigger") for ev in report["isolated"])
    if strongest and not report["isolated"]:
        events.append(strongest)
    with fits.open(inputs/f"{source}-s{sector}.fits",memmap=True) as ff:
        cube=np.asarray(ff[1].data["FLUX"][row],float)
    yy,xx=np.indices(ap.shape)
    edge=(yy==0)|(xx==0)|(yy==ap.shape[0]-1)|(xx==ap.shape[1]-1)
    sky=np.nanmedian(cube[:,edge],axis=1)
    sub=cube-sky[:,None,None]
    tx,ty=rec["target_xy"]
    core=(xx==int(round(tx)))&(yy==int(round(ty)))
    square=(abs(xx-round(tx))<=1)&(abs(yy-round(ty))<=1)
    local=(xx-tx)**2+(yy-ty)**2<=2.5**2
    metrics=[]
    for ev in events:
        in_,left,right=event_masks(t,ev["epoch_btjd"],ev["duration_min"]/1440,ev.get("period_days"))
        row_out=dict(ev,n_in=int(in_.sum()),n_left=int(left.sum()),n_right=int(right.sum()))
        if not (in_.any() and left.any() and right.any()):
            metrics.append(dict(row_out,status="insufficient flanks"))
            continue
        def means(a):
            return np.mean(a[left|right],axis=0),np.mean(a[in_],axis=0)
        # The comparison with the six-hour normalization catches high-flux
        # shoulders masquerading as deep dips when central flux is ordinary.
        row_out["flattened_median_fraction"]={key:float(np.median(resid[mask])) for key,mask in
                                              (("left",left),("in",in_),("right",right))}
        row_out["flattened_mean_fraction"]={key:float(np.mean(resid[mask])) for key,mask in
                                            (("left",left),("in",in_),("right",right))}
        shoulder,center=means(sub)
        deficit=shoulder-center
        row_out["depth_fraction"]={}
        for name,mask in (("central_pixel",core),("aperture",ap),("three_by_three",square)):
            row_out["depth_fraction"][name]=float(np.sum(deficit[mask])/np.sum(shoulder[mask]))
        bg_out,bg_in=means(sky)
        xout,xin=means(cx);yout,yin=means(cy)
        peak=np.nanargmax(np.where(local,deficit,-np.inf))
        py,px=np.unravel_index(peak,ap.shape)
        row_out.update(status="measured",sky_change_e_s_per_pixel=float(bg_in-bg_out),
                       centroid_delta_x_pix=float(xin-xout),centroid_delta_y_pix=float(yin-yout),
                       strongest_local_deficit_pixel_xy=[int(px),int(py)],
                       strongest_pixel_distance_to_WD_pix=float(np.hypot(px-tx,py-ty)),
                       difference_image_e_s=np.round(deficit,3).tolist())
        metrics.append(row_out)
    return dict(tic=tic,gaia=source,sector=sector,target_xy=rec["target_xy"],
                n_periodic=len(report["periodic"]),n_isolated=len(report["isolated"]),image_checks=metrics)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("inputs","curves","search_code","out"):
        p.add_argument(name,type=Path)
    args=p.parse_args()
    search=module(args.search_code)
    obs=json.loads((args.curves/"ffi_measurements.json").read_text())
    rows=[describe(q,search,args.inputs,args.curves) for q in obs if q["status"]=="included"]
    args.out.write_text(json.dumps(rows,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
