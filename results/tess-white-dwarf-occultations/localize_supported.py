"""Fit the WD image and a local sky plane to the mask-supported isolated boxes on SPOC pixels."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from ffi_acquire import propagated
from ffi_images import event_masks
from replication import direct
from short_pixels import acquire
from vet_extension_pixels import sha
from wd_search import mad


def measure(event, star, filename, cache):
    path = cache / filename
    url = direct(filename)
    digest = acquire(path, url)
    epoch = event["epoch_btjd"]
    width = event["duration_min"]/1440
    sector = int(filename.split("-")[1][1:])
    with fits.open(cache/f"{star['GaiaEDR3']}-s{sector}.fits") as h:
        t, corr = (h[1].data[key].astype(float) for key in ("TIME","TIMECORR"))
        valid = np.isfinite(t+corr)
        sc_ffi, correction = (t-corr)[valid], corr[valid]
    with fits.open(path) as h:
        if int(h[0].header["TICID"])!=int(star["tic_requested"]) or int(h[0].header["SECTOR"])!=sector:
            raise ValueError("Target pixel identity mismatch")
        if abs(h[1].header["TIMEDEL"]*86400-120)>.1:raise ValueError("Unexpected cadence")
        d = h[1].data
        time_sc = d["TIME"].astype(float)-d["TIMECORR"]
        t = time_sc+np.interp(time_sc,sc_ffi,correction)
        valid = (d["QUALITY"]==0) & np.isfinite(t)
        t, cube = t[valid], d["FLUX"][valid].astype(float)
        cadence_numbers = d["CADENCENO"][valid].astype(int)
        xy = WCS(h[2].header,naxis=2).world_to_pixel_values(*propagated(star,2000+(epoch+2457000-2451545)/365.25))
    y, x = np.indices(cube.shape[1:])
    # The target's measured shape is frozen from out-of-window stamps, using every 20th exposure.
    oot = np.flatnonzero(abs(t-epoch)>2*width)[::20]
    template = np.nanmedian(cube[oot],axis=0)
    u, v = x-xy[0], y-xy[1]
    disk = (u*u+v*v<=2.5**2) & np.isfinite(template)
    mat = np.stack((template[disk], np.ones(sum(disk.ravel())),u[disk],v[disk]),axis=1)
    weights = np.linalg.pinv(mat)
    levels = cube[:,disk] @ weights.T
    mid, left, right = event_masks(t,epoch,width,None)
    def contrast(series, a, b, c):
        return .5*(series[b].mean(axis=0)+series[c].mean(axis=0))-series[a].mean(axis=0)
    observed = contrast(levels, mid, left, right)
    image = contrast(cube, mid, left, right)
    # Translate the three actual sets of valid cadence numbers together. Keep disjoint
    # off-time triplets with precisely the same retained counts, gaps and spacecraft phase.
    counts = [mask.sum() for mask in (mid,left,right)]
    target_indices = [cadence_numbers[mask] for mask in (mid,left,right)]
    cadence_step = int(np.ceil(4*width*86400/120))+1
    origin = cadence_numbers[mid][0]
    controls = []
    for n in range((cadence_numbers[0]-origin)//cadence_step,(cadence_numbers[-1]-origin)//cadence_step+1):
        if n==0:continue
        trial = [v+n*cadence_step for v in target_indices]
        samples = [np.minimum(np.searchsorted(cadence_numbers,v),len(cadence_numbers)-1) for v in trial]
        if all(np.array_equal(cadence_numbers[i],v) for i,v in zip(samples,trial)):
            controls.append(contrast(levels,*samples)[0])
    controls = np.asarray(controls)
    core_mask = u*u+v*v<=1.5**2
    return dict(tic=star["tic_requested"],sector=sector,epoch_btjd=epoch,duration_min=width*1440,
                target_xy=list(map(float,xy)),source_url=url,sha256=digest,
                counts_center_left_right=[int(n) for n in counts],pixel_count=int(disk.sum()),
                delta_template_fraction=float(observed[0]),sky_offset_e_s_pixel=float(observed[1]),
                core_aperture_depth=float(image[core_mask].sum() / template[core_mask].sum()),
                controls=len(controls),control_median_fraction=float(np.median(controls)),
                control_mad_sigma_fraction=float(mad(controls)),control_percentiles_2p5_97p5=np.percentile(controls,[2.5,97.5]).tolist(),
                local_template_e_s=np.round(template,4).tolist(),difference_e_s=np.round(image,4).tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("fixed","parent","inventory","cache","output"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    stars={r["tic_requested"]:r for r in csv.DictReader(a.parent.open())}
    roster=list(csv.DictReader(a.inventory.open()))
    out=[]
    for s in json.load(gzip.open(a.fixed))["stars"]:
        for e in s["results"]:
            if e["kind"]=="isolated" and e["checks"]["ffi_clean"]["passes"]:
                row,=[r for r in roster if r["tic"]==s["tic"] and int(r["sector"])==s["sector"] and r["cadence_s"]=="120"]
                out.append(measure(e,stars[s["tic"]],row["filename"],a.cache))
    a.output.write_text(json.dumps(dict(fixed_file_sha256=sha(a.fixed),results=out),indent=1,allow_nan=False)+"\n")
