#!/usr/bin/env python3
"""Measure brightness symmetry and cross-section centroids of profile-gate segments."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from wake_image import ImagePlane


def measure(image, valid, row):
    dx,dy = row["x1"]-row["x0"],row["y1"]-row["y0"]
    extent = np.hypot(dx,dy)
    u = np.array([dx,dy])/extent
    v = np.array([-u[1],u[0]])
    width = row["fwhm_px"]
    across = np.arange(-4*width,4*width+.25,.25)
    nbin = 8
    n = max(4,round(extent/nbin))
    along = (np.arange(n*nbin)+.5)/(n*nbin)*extent
    x = row["x0"]+u[0]*along+v[0]*across[:,None]
    y = row["y0"]+u[1]*along+v[1]*across[:,None]
    if not np.all(map_coordinates(valid,[y,x],order=0,mode="constant",cval=0)):
        return dict(usable=0)
    dat = map_coordinates(image,[y,x],order=1,mode="constant",cval=0)
    # Sum cross-window surface brightness along each bin. Subtract the median of
    # the two side bands for a robust local zero point independently in each bin.
    binned = dat.reshape(len(across),nbin,n).mean(axis=2)
    binned -= np.median(binned[np.abs(across)>=2.5*width],axis=0)
    center = np.abs(across)<=.5*width
    flux = binned[center].mean(axis=0)
    pos = np.maximum(flux,0.)
    if pos.sum()<=0:
        return dict(usable=0)
    coords = (np.arange(nbin)+.5)/nbin
    cross = np.abs(across)<=1.5*width
    weights = np.maximum(binned[cross],0)
    norm = weights.sum(axis=0)
    good = (norm>0) & (flux>0)
    cent = np.full(nbin,np.nan)
    cent[good] = (across[cross,None]*weights).sum(axis=0)[good]/norm[good]
    smoothness = None
    if good.sum()>=5:
        coef = np.polyfit(coords[good],cent[good],deg=1)
        smoothness = round(float(np.sqrt(np.mean((cent[good]-np.polyval(coef,coords[good]))**2))
                                 /width),4)
    return dict(usable=1,positive_flux_centroid_0_to_1=round(float(pos@coords/pos.sum()),4),
                brightest_bin_0_to_7=int(np.argmax(flux)),
                inner_half_fraction=round(float(pos[2:6].sum()/pos.sum()),4),
                reversal_asymmetry=round(float(np.abs(pos-pos[::-1]).sum()/(2*pos.sum())),4),
                centerline_scatter_in_fwhm=smoothness,
                mean_native_by_bin=" ".join(f"{z:.7g}" for z in flux),
                transverse_centroid_px=" ".join(f"{z:.3f}" if np.isfinite(z) else "NA" for z in cent))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--runs",type=Path,required=True)
    p.add_argument("--images",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    data=[]
    for path in sorted(a.runs.glob("*.json")):
        if path.name.split("_")[-1] not in ("invert.json","search.json"):
            continue
        r=json.loads(path.read_text())
        rows=[(i,d) for i,d in enumerate(r["detections"]) if d["profile_gate"]]
        if not rows: continue
        tile=path.name.split("_")[0]
        sign = "negative" if r["params"]["invert"] else "positive"
        with ImagePlane(a.images/f"{tile}_VIS.fits") as sci, ImagePlane(a.images/f"{tile}_rms.fits") as wt:
            image = sci.read(0,0,sci.width,sci.height)
            mask = wt.read(0,0,wt.width,wt.height)
        valid = (np.isfinite(image) & (image!=0) & np.isfinite(mask) & (mask>0) & (mask<1e8)).astype("u1")
        if sign=="negative": image = -image
        for i,s in rows:
            data.append(dict(id=f"{tile}_{sign}_{i}",tile=tile,sign=sign,**measure(image,valid,s)))
    cols=list(dict.fromkeys(k for r in data for k in r))
    with a.out.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=cols,lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


if __name__=="__main__":
    main()
