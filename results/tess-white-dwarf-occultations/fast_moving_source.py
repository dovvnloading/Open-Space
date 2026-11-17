"""Measure the moving residual in the fast TPF that gave the largest inferred central loss."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.special import ndtr

from pixel_evidence import profile_center
from short_pixels import acquire


def profile(x,y,cx,cy,sigma):
    return ((ndtr((x+.5-cx)/sigma)-ndtr((x-.5-cx)/sigma)) *
            (ndtr((y+.5-cy)/sigma)-ndtr((y-.5-cy)/sigma)))


def run(record, trigger, cache):
    epoch=trigger["epoch_btjd"]
    width=trigger["duration_min"]/1440
    path=cache/record["url"].split("/")[-1]
    if acquire(path,record["url"])!=record["sha256"]:raise ValueError("Fast TPF hash mismatch")
    with fits.open(path) as h:
        d=h[1].data
        ok=(d["QUALITY"]==0)&np.isfinite(d["TIME"])
        t=d["TIME"][ok].astype(float)
        cube=d["FLUX"][ok].astype(float)
    # The exposure intervals lie on both sides of the moving object's passage through the stamp.
    reference_mask=(abs(t-epoch)>.28)&(abs(t-epoch)<.55)
    reference=np.median(cube[reference_mask],axis=0)
    bins=[]
    for offset in np.arange(-.10,.21,.01):
        use=abs(t-epoch-offset)<.005
        if not use.any():continue
        delta=np.median(cube[use],axis=0)-reference
        xy=profile_center(delta)
        if xy is not None:
            bins.append(dict(btjd=float(np.median(t[use])),x=float(xy[0]),y=float(xy[1]),
                             peak_e_s=float(np.max(delta-np.median(delta))),points=int(use.sum())))
    # Only bright peaks away from all stamp edges constrain the track.
    use=[r for r in bins if r["peak_e_s"]>15 and min(r["x"],r["y"],10-r["x"],10-r["y"])>1]
    pos=np.array([[r["x"],r["y"]] for r in use])
    design=np.array([[1,r["btjd"]-epoch] for r in use])
    coefficients=np.linalg.lstsq(design,pos,rcond=None)[0]
    rms=np.sqrt(np.mean(np.sum((design@coefficients-pos)**2,axis=1)))
    x0,y0=record["target_xy"]
    y,x=np.indices(reference.shape)
    disk=np.hypot(x-x0,y-y0)<=2.5
    ap=np.hypot(x-x0,y-y0)<=1.5
    center=abs(t-epoch)<=.5*width
    image=np.mean(cube[center],axis=0)-reference
    xy=np.stack((np.ones(center.sum()),t[center]-epoch),axis=1)@coefficients
    # Fit the time-integrated moving template and the stationary empirical stellar
    # image in the same local disk; free the sky, its two slopes, and both amplitudes.
    models=[]
    for sigma in (.65,.85,1.05):
        moving=np.mean([profile(x,y,*point,sigma) for point in xy],axis=0)
        mat=np.stack((reference[disk],moving[disk],np.ones(disk.sum()),x[disk]-x0,y[disk]-y0),axis=1)
        fit=np.linalg.lstsq(mat,image[disk],rcond=None)[0]
        models.append(dict(sigma_pix=sigma,wd_loss_vs_quiet=float(-fit[0]),moving_flux_e_s=float(fit[1])))
    return dict(tic=record["tic"],sector=record["sector"],source_url=record["url"],sha256=record["sha256"],
                trigger=trigger,target_xy=[x0,y0],reference_time_offset_days=[.28,.55],
                quiet_samples=int(reference_mask.sum()),central_samples=int(center.sum()),
                centers=bins,track_centers=len(use),track_xy_at_epoch=coefficients[0].tolist(),
                vx_vy_pix_per_day=coefficients[1].tolist(),track_rms_pix=float(rms),
                aperture_excess_vs_quiet=float(image[ap].sum()/reference[ap].sum()),moving_plus_wd_fits=models)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("fast",type=Path);p.add_argument("cache",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args()
    rows=json.load(gzip.open(a.fast))
    r,e=max(((r,e) for r in rows for e in r["events"] if e["status"]=="measured"),
            key=lambda pair:pair[1]["wd_image"]["center_loss_vs_context"])
    a.out.write_text(json.dumps(run(r,e["trigger"],a.cache),indent=1,allow_nan=False)+"\n")
