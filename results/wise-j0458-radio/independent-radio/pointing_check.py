"""Evaluate matched Stokes I/V apertures in native SURF per-pointing images."""
import csv
import hashlib
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from verify import ROOT, fit_at, position


def main():
    with (ROOT/"single-pointing-products.csv").open() as f:
        entries=list(csv.DictReader(f))
    rows=[]
    for r in entries:
        path=ROOT/r["local_product"]
        if hashlib.sha256(path.read_bytes()).hexdigest()!=r["crop_sha256"]:
            raise ValueError("Product differs: "+path.name)
        with fits.open(path) as hdus:
            h=hdus[0].header
            data=hdus[0].data.squeeze().astype(float)*1000
        if h["CRVAL3"]!=float(r["stokes_fits_code"]) or h["BUNIT"].lower()!="jy/beam":
            raise ValueError("Unexpected Stokes plane or units")
        w=WCS(h).celestial
        x,y=w.world_to_pixel_values(*position(Time(h["DATE-OBS"]).mjd))
        jac=w.pixel_scale_matrix*3600
        steps=[(i,j) for i in range(-23,24,6) for j in range(-23,24,6)
               if np.hypot(i,j)>8]
        b=(6.,6.)
        controls=np.array([fit_at(data,x+i,y+j,jac,b) for i,j in steps])
        median=np.median(controls)
        mad=1.4826*np.median(np.abs(controls-median))
        a=fit_at(data,x,y,jac,b)
        rows.append(dict(local_product=path.name,date=h["DATE-OBS"],forced_mjy=a,
                         num_control=len(controls),control_median_mjy=median,
                         control_mad_scale_mjy=mad,control_std_mjy=controls.std(ddof=1),
                         number_control_lower=int(np.sum(controls<a))))
    with (ROOT/"single-pointing-measurements.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__=="__main__":
    main()
