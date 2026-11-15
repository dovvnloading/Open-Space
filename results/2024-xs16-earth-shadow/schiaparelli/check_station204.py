"""Replay the pre-ingress display-image measurement, downloading the six ESA JPEGs.

From this directory run: python check_station204.py
Requires numpy, scipy and Pillow. Downloaded files are verified against the manifest.
"""
import hashlib
import io
import json
import math
import urllib.request
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates
from scipy.optimize import root


def project(h, ra, dec):
    """Gnomonic projection followed by inversion of the FITS forward SIP distortion."""
    longitude = np.deg2rad(ra-h["CRVAL1"])
    d, d0 = np.deg2rad([dec,h["CRVAL2"]])
    cosc = np.sin(d)*np.sin(d0)+np.cos(d)*np.cos(d0)*np.cos(longitude)
    plane = np.rad2deg([np.cos(d)*np.sin(longitude)/cosc,
                       (np.sin(d)*np.cos(d0)-np.cos(d)*np.sin(d0)*np.cos(longitude))/cosc])
    cd = [[h["CD1_1"],h["CD1_2"]],[h["CD2_1"],h["CD2_2"]]]
    uv = np.linalg.solve(cd,plane)
    def distorted(xy):
        result = xy.copy()
        for axis,prefix in enumerate(("A","B")):
            for key,value in h.items():
                if key.startswith(prefix+"_") and key[-1].isdigit():
                    _,i,j = key.split("_")
                    result[axis] += value*xy[0]**int(i)*xy[1]**int(j)
        return result-uv
    x,y = root(distorted,uv).x+[h["CRPIX1"]-1,h["CRPIX2"]-1]
    return x,h["NAXIS2"]-1-y  # published display starts at the top


def score(a, x, y):
    """A 3-pixel Gaussian-weighted aperture minus a local 7-12 pixel sky annulus."""
    dx,dy = np.meshgrid(np.arange(-12,13),np.arange(-12,13))
    radius = np.hypot(dx,dy)
    cut = map_coordinates(a,[y+dy,x+dx],order=1)
    aperture = radius<=3
    weights = np.exp(-radius[aperture]**2/(2*1.5**2))
    return float(np.average(cut[aperture],weights=weights)-np.median(cut[(radius>=7)&(radius<=12)]))


def calculate(directory=None):
    base = Path(__file__).resolve().parent
    frames = json.loads((base/"station204_manifest.json").read_text())
    samples = []; per_frame = []
    offsets = [(ix,iy) for ix in range(-100,101,10) for iy in range(-100,101,10) if 24<=math.hypot(ix,iy)<=100]
    for frame in frames:
        if directory is None:
            with urllib.request.urlopen(frame["url"],timeout=90) as response:
                data=response.read()
        else:
            data=(Path(directory)/frame["local_name"]).read_bytes()
        if hashlib.sha256(data).hexdigest()!=frame["sha256"]:
            raise ValueError("Image checksum changed: "+frame["filename"])
        display=np.asarray(Image.open(io.BytesIO(data)).convert("L"),dtype=float)
        x,y=project(frame["wcs"],frame["ra_horizons"],frame["dec_horizons"])
        # Coordinates and translated tracks were fixed from the ephemeris before inspecting the stack.
        per_frame.append(score(display,x,y))
        samples.append([score(display,x+ix,y+iy) for ix,iy in offsets])
    controls=np.mean(samples,axis=0)
    return {"six_fixed_position_scores":per_frame,"mean_on_track_display_levels":np.mean(per_frame),
            "translated_track_count":len(controls),"translated_tracks_equal_or_above":int(sum(controls>=np.mean(per_frame))),
            "translated_track_median":np.median(controls),"translated_track_max":np.max(controls)}


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--image-dir",help="Use already-downloaded frame1.jpg through frame6.jpg")
    args=parser.parse_args()
    print(json.dumps(calculate(args.image_dir),indent=2))
