#!/usr/bin/env python3
"""Screen one science stamp or measure incremental morphology recovery on real backgrounds."""
import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from wake_filter import FWHM_SIGMA, axes, compatible, detect
from wake_image import ImagePlane


def digest(path):
    with open(path, "rb") as stream:
        return sha256(stream.read()).hexdigest()


def values(s):
    v = [float(x) for x in s.split(",")]
    if not all(np.isfinite(x) and x > 0 for x in v):
        raise argparse.ArgumentTypeError("expected positive finite comma-separated numbers")
    return v


def add_wake(img, truth, intensity):
    """Observed Gaussian transverse FWHM, uniform axial brightness, square-pixel integration."""
    yy, xx = np.indices(img.shape, dtype="f4")
    dx, dy = xx-truth["x"], yy-truth["y"]
    u, v = axes(np.radians(truth["pa_image_deg"]))
    signal = np.zeros_like(img)
    width = truth["fwhm_px"]
    # Normalize the continuous model's mean intensity in its transverse FWHM strip.
    mean_to_peak = np.sqrt(np.pi/(4*np.log(2))) * math.erf(np.sqrt(np.log(2)))
    for sub_x in [-.375,-.125,.125,.375]:
        for sub_y in [-.375,-.125,.125,.375]:
            a = (dx+sub_x)*u[0]+(dy+sub_y)*u[1]
            c = (dx+sub_x)*v[0]+(dy+sub_y)*v[1]
            signal += (np.abs(a) <= truth["length_px"]/2)*np.exp(-.5*(c/(width/FWHM_SIGMA))**2)
    return img + signal*(intensity/(16*mean_to_peak))


def offset(row, x0, y0):
    for key in ("x", "x0", "x1"):
        row[key] += x0
    for key in ("y", "y0", "y1"):
        row[key] += y0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["search", "inject"])
    p.add_argument("image", type=Path)
    p.add_argument("--hdu", type=int, default=0)
    p.add_argument("--scale", type=float, required=True, help="arcsec per square pixel")
    p.add_argument("--window", type=int, nargs=4, metavar=("X", "Y", "W", "H"))
    p.add_argument("--mask", type=Path, help="aligned science mask or weight image")
    p.add_argument("--mask-hdu", type=int, default=0)
    p.add_argument("--mask-type", choices=["bad", "weight", "rms"], default="bad")
    p.add_argument("--ignore-zero", action="store_true", help="science zero is a no-data sentinel")
    p.add_argument("--lengths", type=values, default=[6,10], help="filter lengths in arcsec")
    p.add_argument("--widths", type=values, default=[.20], help="observed filter FWHM in arcsec")
    p.add_argument("--threshold", type=float, default=6., help="empirical per-template score")
    p.add_argument("--min-length", type=float, default=3., help="minimum measured extent, arcsec")
    p.add_argument("--invert", action="store_true", help="negative-image artifact control")
    p.add_argument("--hosts", type=Path, help="CSV id,x,y of galaxies, stamp pixels, zero-based")
    p.add_argument("--host-radius", type=float, default=1., help="end-to-centroid arcsec")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--mu", type=values, default=[24,25,26,27], help="mean AB mag/arcsec² in FWHM strip")
    p.add_argument("--trials", type=int, default=4, help="attempts per brightness-length-width cell")
    p.add_argument("--seed", type=int, default=12023)
    p.add_argument("--patch", type=int, default=192, help="injection trial square side in pixels")
    cal = p.add_mutually_exclusive_group()
    cal.add_argument("--zeropoint", type=float, help="point-source AB ZP for one native pixel unit")
    cal.add_argument("--mjy-sr", action="store_true", help="native pixels in MJy/sr")
    cal.add_argument("--njy", action="store_true", help="native pixels in nJy/pixel")
    args = p.parse_args()
    if args.scale <= 0 or args.trials < 1 or args.min_length <= 0:
        p.error("scale, lengths and trial count must be positive")
    if args.mode == "inject" and (args.invert or not (
            args.zeropoint is not None or args.mjy_sr or args.njy)):
        p.error("injection requires a calibration and a positive science image")
    science = ImagePlane(args.image, args.hdu)
    x0, y0, w, h = args.window or (0, 0, science.width, science.height)
    image = science.read(x0, y0, w, h)
    dims = (science.width, science.height)
    science.close()
    good = np.isfinite(image)
    if args.ignore_zero:
        good &= image != 0
    inputs = {"science": dict(name=args.image.name, sha256=digest(args.image), hdu=args.hdu)}
    if args.mask:
        reader = ImagePlane(args.mask, args.mask_hdu)
        if (reader.width, reader.height) != dims:
            raise ValueError("mask dimensions differ from science image")
        mask = reader.read(x0, y0, w, h)
        reader.close()
        good &= np.isfinite(mask)
        if args.mask_type == "bad":
            good &= mask == 0
        else:
            good &= mask > 0
            if args.mask_type == "rms":
                good &= mask < 1e8
        inputs["mask"] = dict(name=args.mask.name, sha256=digest(args.mask), hdu=args.mask_hdu)
    lengths = np.array(args.lengths) / args.scale
    widths = np.array(args.widths) / args.scale
    def screen(data, valid):
        return detect(data, valid, lengths, widths, args.scale, args.threshold, args.min_length)
    result = dict(mode=args.mode, inputs=inputs, params={k:str(v) if isinstance(v,Path) else v
                  for k,v in vars(args).items() if k not in ("out","image","mask","hosts")},
                  window_xywh=[x0,y0,w,h], valid_pixels=int(good.sum()),
                  valid_geometric_deg2=float(good.sum()*(args.scale/3600)**2))
    if args.mode == "search":
        hosts = []
        if args.hosts:
            inputs["host_catalog"] = dict(name=args.hosts.name, sha256=digest(args.hosts))
            with args.hosts.open() as stream:
                hosts = list(csv.DictReader(stream))
        detections = screen(-image if args.invert else image, good)
        for row in detections:
            offset(row, x0, y0)
            if args.hosts:
                row["endpoint_galaxies"] = []
                for host in hosts:
                    distance = min(math.hypot(row["x"+str(i)]-float(host["x"]),
                                              row["y"+str(i)]-float(host["y"])) for i in (0,1))
                    if distance*args.scale <= args.host_radius:
                        row["endpoint_galaxies"].append(dict(id=host["id"], distance_arcsec=distance*args.scale))
        result["detections"] = detections
    else:
        rng = np.random.default_rng(args.seed)
        trials = []
        size = min(args.patch, h, w)
        margin = math.ceil(max(lengths)*.60+4*max(widths))
        if size <= 2*margin+2:
            raise ValueError("trial patch too small for template and margins")
        for length in lengths:
            for width in widths:
                for index in range(args.trials):
                    # Paired magnitudes share sites and orientations. Different sites need
                    # not be independent if trial patches overlap. Draw from the whole stamp.
                    px = int(rng.integers(w-size+1))
                    py = int(rng.integers(h-size+1))
                    region = np.s_[py:py+size,px:px+size]
                    src, valid = image[region], good[region]
                    truth = dict(x=float(rng.uniform(margin,size-1-margin)), y=float(
                        rng.uniform(margin,size-1-margin)), pa_image_deg=float(rng.uniform(0,180)),
                        fwhm_px=float(width), length_px=float(length))
                    before = screen(src,valid)
                    pre = any(compatible(truth,d) for d in before)
                    for mu in args.mu:
                        if args.mjy_sr:
                            intensity = 3631*10**(-.4*mu)*(180*3600/np.pi)**2/1e6
                        else:
                            zp = 2.5*np.log10(3631e9) if args.njy else args.zeropoint
                            intensity = args.scale**2 * 10**(-.4*(mu-zp))
                        output = screen(add_wake(src,truth,intensity),valid)
                        match = [d for d in output if compatible(truth,d)]
                        trials.append(dict(index=index, x=x0+px+truth["x"], y=y0+py+truth["y"],
                            patch_xy=[x0+px,y0+py], angle_deg=truth["pa_image_deg"],
                            length_arcsec=float(length*args.scale), fwhm_arcsec=float(width*args.scale),
                            mu_ab_arcsec2=mu, mean_native_per_pixel=float(intensity),
                            preexisting_match=pre, detected=bool(match), new_match=bool(match) and not pre,
                            matched_score=max((d["score"] for d in match),default=None)))
        result["trials"] = trials
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(args.out.suffix+".tmp")
    tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    tmp.replace(args.out)


if __name__ == "__main__":
    main()
