"""Compare both predicted tracks in selected HST exposures at JWST threshold peaks.

The search around each predicted location allows 0.15 arcsec. The significance
values are descriptive local MAD scores, with no correction for searching offsets.
"""
import argparse
import csv
import json
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy import ndimage

from fetch_hst_images import HERE, main as check_inputs
BASE = HERE.parent
NIR = BASE / "nircam-hd53143-hd38206/outputs"
DISKS = BASE / "hd21997-hd202628"
sys.path.insert(0, str(BASE.parent / "jwst-nearby-companions/vetting"))
from other_tracks import offset, parpos

HOSTS = {"HD53143": "5479222240596469632", "HD38206": "2967287596507432704",
         "HD21997": "5084963312478488960", "HD202628": "6577318075351824512"}


def records(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def peaks(field):
    peaks = []
    if field in ("HD53143", "HD38206"):
        paths = NIR.glob(field.lower() + "*_peaks.csv")
    else:
        paths = [DISKS / "measurements" / (field+"_peaks.csv")]
        for row in records(DISKS / "measurements/auxiliary_peaks.csv"):
            if row["field"] == field:
                peaks.append((row["file"], row["id"], float(row["mjd"]), float(row["east"]), float(row["north"])))
    for path in sorted(paths):
        meta = json.loads(path.with_name(path.name.replace("_peaks.csv", "_measurement.json")).read_text())
        for i, row in enumerate(records(path), 1):
            peaks.append((path.name, i, meta["expstart_mjd"], float(row["east_arcsec"]), float(row["north_arcsec"])))
    return peaks


def load_frame(path, star):
    with fits.open(path) as hd:
        primary = hd[0].header
        if "SNR_MAP" in hd:
            hdr = hd["DATA_INFORMATION"].data[0]
            img = hd["SNR_MAP"].data[0].astype(float)
            sens = hd["SENSITIVITY_MAP"].data[0]
            img[(sens <= 0) | ~np.isfinite(sens)] = np.nan
            scale = primary["PIXSCALE"]
            return dict(name=path.name, native=img, time=float(hdr["EXPSTART"]) if "EXPSTART" in hdr.array.names else float(hdr[5]),
                        mat=np.diag([-scale, scale]), center=np.array([(img.shape[1]-1)/2,(img.shape[0]-1)/2]),
                        method="ALICE north-combined SNR", instrument="NICMOS"+primary["FILTER"])
        sh = hd["SCI"].header
        wcs = WCS(sh, fix=False)
        mjd = primary.get("EXPSTART", primary.get("TEXPSTRT"))
        image = hd["SCI"].data.astype(float)
        stis = hd[0].header["INSTRUME"] == "STIS"
        if stis:
            center = np.array([sh["CRPIX1"]-1, sh["CRPIX2"]-1])
            valid = hd["DQ"].data == 0
            image /= primary["TEXPTIME"]
        else:
            center = np.array(wcs.world_to_pixel_values(*parpos(star, Time(mjd, format="mjd"))))
            valid = hd["WHT"].data > 0.25 * np.nanmax(hd["WHT"].data)
        valid &= np.isfinite(image)
        image[~valid] = np.nan
        return dict(name=path.name, native=image, time=mjd, mat=wcs.pixel_scale_matrix*3600,
                    center=center, method="STIS detector median subtraction" if stis else "ACS high-pass",
                    instrument="STIS" if stis else "ACS")


def filter_frames(frames):
    for f in frames:
        raw = f["native"]
        if f["instrument"].startswith("NICMOS"):
            f["score"] = raw
            continue
        valid = np.isfinite(raw)
        if f["instrument"] == "STIS":
            # These rolls share the aperture and detector subarray. Use independent
            # orientations so a point at the tested separations moves off itself.
            a = np.rad2deg(np.arctan2(f["mat"][0,1], f["mat"][1,1]))
            others = []
            for k in frames:
                if k["instrument"] != "STIS" or k is f:
                    continue
                b = np.rad2deg(np.arctan2(k["mat"][0,1], k["mat"][1,1]))
                if abs((a-b+180)%360-180) > 8:
                    others.append(k["native"])
            if not others:
                raise ValueError("No independent HST roll")
            ref = np.nanmedian(others, axis=0)
            raw = raw - ref
            valid &= np.isfinite(ref)
        z = np.where(valid,raw,0.)
        # Compact positive responses on scales of a few native pixels.
        residual = ndimage.gaussian_filter(z,1) - ndimage.gaussian_filter(z,4)
        residual[~ndimage.minimum_filter(valid,5)] = np.nan
        f["score"] = residual


def measure(frame, xy):
    x,y = frame["center"] + np.linalg.solve(frame["mat"], xy)
    scale = np.sqrt(abs(np.linalg.det(frame["mat"])))
    image = frame["score"]
    m = int(np.ceil(.48/scale))+2
    i,j = int(round(x)),int(round(y))
    if i-m<0 or j-m<0 or i+m>=image.shape[1] or j+m>=image.shape[0]:
        return None
    small = image[j-m:j+m+1, i-m:i+m+1]
    yy,xx = np.mgrid[j-m:j+m+1,i-m:i+m+1]
    rad = np.hypot(xx-x,yy-y)*scale
    signal = (rad <= .15) & np.isfinite(small)
    ring = small[(rad >= .22) & (rad < .48) & np.isfinite(small)]
    if signal.sum() < 3 or len(ring) < 18:
        return None
    med = np.median(ring)
    sigma = 1.4826*np.median(abs(ring-med))
    if not sigma>0:
        return None
    if frame["instrument"].startswith("NICMOS"):
        med,sigma=0.,1.
    py,px = np.unravel_index(np.argmax(np.where(signal,small,-np.inf)),small.shape)
    de,dn = frame["mat"] @ np.array([xx[py,px]-x,yy[py,px]-y])
    return dict(forced_score=(image[j,i]-med)/sigma if np.isfinite(image[j,i]) else "",
                max_score=(small[py,px]-med)/sigma, peak_offset_e=de,peak_offset_n=dn)


def main(cache):
    inputs = check_inputs(cache)
    out = []
    for field, source in HOSTS.items():
        cat = (NIR / (field+"-gaia.csv")) if field in ("HD38206", "HD53143") else (DISKS / "catalogs" / (field+"_GaiaDR3.csv"))
        star = next(row for row in records(cat) if row["source_id"] == source)
        frames = [load_frame(cache/r["name"], star) for r in inputs if r["field"]==field]
        filter_frames(frames)
        for f in frames:
            pos = parpos(star,Time(f["time"],format="mjd"))
            for product,peak,mjd,e,n in peaks(field):
                shift = offset(parpos(star,Time(mjd,format="mjd")),pos)
                for track, xy in [("comoving",np.array([e,n])), ("fixed",np.array([e,n])+shift)]:
                    result = measure(f,xy)
                    out.append(dict(field=field,jwst_product=product,peak=peak,hst_product=f["name"],mjd=f["time"],
                                    method=f["method"],track=track,expected_e=xy[0],expected_n=xy[1],covered=int(result is not None),
                                    **(result or dict(forced_score="",max_score="",peak_offset_e="",peak_offset_n=""))))
    with (HERE/"hst_tracks.csv").open("w") as stream:
        writer=csv.DictWriter(stream,out[0]);writer.writeheader();writer.writerows(out)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache",type=Path)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        main(p.parse_args().cache)
