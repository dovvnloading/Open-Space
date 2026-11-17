#!/usr/bin/env python3
"""Combine per-stamp receipts and profiles, preserving matched rows of both signs."""
import argparse
import csv
import gzip
import io
import json
import math
from pathlib import Path
from collections import Counter, defaultdict

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from scipy.ndimage import distance_transform_edt, map_coordinates
from scipy.spatial import cKDTree

from wake_image import ImagePlane


def save_rows(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    s = io.StringIO()
    wr = csv.DictWriter(s,fieldnames=fields,lineterminator="\n")
    wr.writeheader()
    wr.writerows(rows)
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(s.getvalue().encode(),mtime=0))
    else:
        path.write_text(s.getvalue())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--runs",type=Path,required=True)
    p.add_argument("--images",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=True,parents=True)
    coverage,segments,trials,short = [],[],[],[]
    centers = []
    for f in sorted(a.runs.glob("*_search.json")):
        tile = f.name.split("_")[0]
        d = json.loads(f.read_text())
        if d["mode"] != "search" or d["params"]["invert"]:
            raise ValueError(f)
        full = a.images / d["inputs"]["science"]["name"]
        h = fits.getheader(full)
        if h["CTYPE1"] != "RA---TAN" or h["CTYPE2"] != "DEC--TAN":
            raise ValueError("Unexpected survey projection")
        wcs = WCS(h,naxis=2)
        _,_,w,height = d["window_xywh"]
        with ImagePlane(full) as vis_reader, ImagePlane(a.images/d["inputs"]["mask"]["name"]) as rms_reader:
            raw = vis_reader.read(0,0,w,height)
            rms = rms_reader.read(0,0,w,height)
        ok = np.isfinite(raw) & (raw != 0) & np.isfinite(rms) & (rms>0) & (rms<1e8)
        if int(ok.sum()) != d["valid_pixels"]:
            raise ValueError("mask pixel-count mismatch")
        sc = d["params"]["scale"]
        halo = math.ceil(.9*max(d["params"]["lengths"])/sc + 4*max(d["params"]["widths"])/sc) + 2
        # Sufficient radius to fit the full measured stripe at every filter angle.
        dist = distance_transform_edt(np.pad(ok,1))[1:-1,1:-1]
        info = dict(tile=tile,vis_sha256=d["inputs"]["science"]["sha256"],
                    rms_sha256=d["inputs"]["mask"]["sha256"], nx=w,ny=height,
                    valid_pixels=int(ok.sum()),all_angles_center_pixels=int((dist>halo).sum()),
                    center_margin_px=halo, square_arcsec_per_pixel=sc*sc)
        center = wcs.pixel_to_world_values((w-1)/2,(height-1)/2)
        centers.append(tuple(float(x) for x in center))
        inv = a.runs / f"{tile}_invert.json"
        runs = [(1,d)]
        if inv.exists():
            runs.append((-1,json.loads(inv.read_text())))
        for sign,run in runs:
            label = "positive" if sign==1 else "negative"
            info[label+"_raw"] = len(run["detections"])
            info[label+"_profile"] = sum(s["profile_gate"] for s in run["detections"])
            for i,s in enumerate(run["detections"]):
                payload = dict(tile=tile,sign=sign,rank=i+1,**s)
                for idx in ("","0","1"):
                    ra,dec = wcs.pixel_to_world_values(s["x"+idx],s["y"+idx])
                    payload["ra"+idx] = round(float(ra),8)
                    payload["dec"+idx] = round(float(dec),8)
                segments.append(payload)
                if sign==1 and s["profile_gate"]:
                    short.append(payload)
        coverage.append(info)
    sample_file = a.runs / "injection_sample.csv"
    with sample_file.open() as stream:
        samples = {r["tile"]:r["selection"] for r in csv.DictReader(stream)}
    for f in sorted(a.runs.glob("*_inject.json")):
        tile = f.name.split("_")[0]
        run = json.loads(f.read_text())
        sc = run["params"]["scale"]
        im_path = a.images / run["inputs"]["science"]["name"]
        rms_path = a.images / run["inputs"]["mask"]["name"]
        with ImagePlane(im_path) as vis_reader, ImagePlane(rms_path) as rms_reader:
            w,h = vis_reader.width,vis_reader.height
            im = vis_reader.read(0,0,w,h)
            rm = rms_reader.read(0,0,w,h)
        valid = (np.isfinite(im) & (im != 0) & np.isfinite(rm) & (rm>0) & (rm<1e8)).astype("u1")
        for r in run["trials"]:
            r = dict(r)
            ang = np.radians(r["angle_deg"])
            length,width = r["length_arcsec"]/sc,r["fwhm_arcsec"]/sc
            along = np.arange(-length/2,length/2+.25,.5)
            cross = np.arange(-4*width,4*width+.25,.5)[:,None]
            x = r["x"]+along*np.cos(ang)-cross*np.sin(ang)
            y = r["y"]+along*np.sin(ang)+cross*np.cos(ang)
            region = map_coordinates(valid,[y,x],order=0,mode="constant",cval=0)
            r["strip_valid"] = bool(region.min())
            r["selection"] = samples[tile]
            r["patch_x"],r["patch_y"] = r.pop("patch_xy")
            trials.append(dict(tile=tile,seed=run["params"]["seed"],**r))
    save_rows(a.out/"coverage.csv",coverage)
    save_rows(a.out/"segments.csv.gz",segments)
    save_rows(a.out/"profile_shortlist.csv",short)
    save_rows(a.out/"injections.csv",trials)
    counters = defaultdict(Counter)
    for t in trials:
        k = (t["length_arcsec"], t["fwhm_arcsec"],t["mu_ab_arcsec2"])
        c = counters[k]
        c["attempts"] += 1
        for pred in ("preexisting_match","detected","new_match","new_profile_match"):
            c[pred] += int(t[pred])
    save_rows(a.out/"recovery.csv",[dict(length_arcsec=k[0],fwhm_arcsec=k[1],
              mu_ab_arcsec2=k[2],**v) for k,v in sorted(counters.items())])
    strata = defaultdict(Counter)
    for t in trials:
        for domain in (["all","valid strip"] if t["strip_valid"] else ["all"]):
            key = (t["selection"],domain,t["length_arcsec"],t["fwhm_arcsec"],t["mu_ab_arcsec2"])
            c = strata[key]
            c["attempts"] += 1
            for pred in ("preexisting_match","detected","new_match","new_profile_match"):
                c[pred] += int(t[pred])
    save_rows(a.out/"recovery_by_sample.csv", [dict(selection=k[0],coverage=k[1],
        length_arcsec=k[2],fwhm_arcsec=k[3],mu_ab_arcsec2=k[4],**v)
        for k,v in sorted(strata.items())])
    min_sep = None
    if len(centers)>1:
        ra,dec = np.radians(centers).T
        xyz = np.column_stack((np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)))
        chord = cKDTree(xyz).query(xyz,k=2)[0][:,1].min()
        min_sep = float(np.degrees(2*np.arcsin(chord/2))*3600)
    largest_diagonal = max((math.hypot(r["nx"],r["ny"])*math.sqrt(r["square_arcsec_per_pixel"])
                            for r in coverage), default=0.)
    if min_sep is not None and min_sep <= largest_diagonal:
        raise ValueError("stamps may overlap; geometric area requires union")
    s = dict(searched_tiles=len(coverage),with_valid_pixels=sum(r["valid_pixels"]>0 for r in coverage),
             min_center_separation_arcsec=min_sep,max_stamp_diagonal_arcsec=largest_diagonal,
             valid_area_deg2=sum(r["valid_pixels"]*r["square_arcsec_per_pixel"] for r in coverage)/3600**2,
             full_stripe_center_area_deg2=sum(r["all_angles_center_pixels"]*r["square_arcsec_per_pixel"]
                                              for r in coverage)/3600**2,
             inverse_tiles=sum("negative_raw" in r for r in coverage),
             **{k:sum(r.get(k,0) for r in coverage) for k in ("positive_raw","positive_profile","negative_raw","negative_profile")},
             injection_attempts=len(trials),injection_tiles=len(set(t["tile"] for t in trials)))
    (a.out/"summary.json").write_text(json.dumps(s,indent=2)+"\n")
    print(json.dumps(s))


if __name__ == "__main__":
    main()
