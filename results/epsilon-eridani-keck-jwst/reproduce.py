"""Check the 2001/2023 epsilon Eridani coincidence using KOA and MAST FITS images."""
import argparse
import csv
import hashlib
from pathlib import Path
import shutil
from urllib.request import urlopen
import warnings

import numpy as np
from astropy.io import fits
from astropy.stats import mad_std
from astropy.wcs import WCS, FITSFixedWarning
from scipy import ndimage
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
# Catalog coordinates serve only as search seeds for the nearby bright reference.
REFERENCE_SKY = np.array([53.22629167, -9.45809829]) + np.array([11.30 / 3551.1, 13.45 / 3600])
EAST_NORTH_SCALE = np.array([3600 * np.cos(np.deg2rad(REFERENCE_SKY[1])), 3600])
KECK_VECTOR = np.array([-3.5, 2.2])  # Macintosh Table 3, source 7 minus source 1


def fetch(cache, manifest):
    cache.mkdir(parents=True, exist_ok=True)
    for row in manifest:
        path = cache / row["filename"]
        if not path.exists():
            partial = path.with_suffix(".part")
            with urlopen(row["url"], timeout=120) as response, partial.open("wb") as output:
                shutil.copyfileobj(response, output)
            partial.replace(path)
        with path.open("rb") as data:
            digest = hashlib.file_digest(data, "sha256").hexdigest()
        if digest != row["sha256"]:
            raise ValueError(f"Product version differs from manifest: {path.name} ({digest})")


def aperture(image, xy, radius, inner, outer):
    """Fractional-pixel circle with local median annulus. Returns native pixel sum."""
    x, y = map(float, xy)
    cx, cy = round(x), round(y)
    half = int(np.ceil(outer + 1))
    cut = image[cy-half:cy+half+1, cx-half:cx+half+1]
    if cut.shape != (2*half+1, 2*half+1):
        return np.nan
    dy, dx = np.mgrid[-half:half+1, -half:half+1]
    dx, dy = dx + cx - x, dy + cy - y
    rr = np.hypot(dx, dy)
    annulus = cut[(inner <= rr) & (rr <= outer)]
    if np.isfinite(annulus).mean() < .95:
        return np.nan
    weights = np.zeros_like(cut, dtype=float)
    subpixels = (np.arange(5) + .5) / 5 - .5
    for a in subpixels:
        for b in subpixels:
            weights += (dx+a)**2 + (dy+b)**2 < radius**2
    use = weights > 0
    if not np.isfinite(cut[use]).all():
        return np.nan
    return np.sum((cut[use]-np.nanmedian(annulus))*weights[use])/25


def keck_centroid(image, xy):
    """Moffat profile over a sloping background in a 19 x 19 pixel patch."""
    x, y = map(float, xy)
    cx, cy = round(x), round(y)
    cut = image[cy-9:cy+10, cx-9:cx+10]
    yy, xx = np.mgrid[-9:10, -9:10]
    rr = np.hypot(xx, yy)
    sky = np.median(cut[rr > 5])
    initial = [x-cx, y-cy, 2.7, 2.6, cut[rr < 3].max()-sky, sky, 0, 0]

    def residual(p):
        dx, dy, width, index, amplitude, base, bx, by = p
        psf = (1 + ((xx-dx)**2 + (yy-dy)**2)/width**2)**(-index)
        return (base + bx*xx + by*yy + amplitude*psf - cut).ravel()

    bounds = ([-4,-4,.9,1.3,0,-np.inf,-np.inf,-np.inf],
              [4,4,8,7,np.inf,np.inf,np.inf,np.inf])
    result = least_squares(residual, initial, bounds=bounds).x
    return np.array([cx, cy]) + result[:2]


def keck_stack(files):
    headers = [fits.getheader(f) for f in files]
    observations = []
    for f in files:
        frame = fits.getdata(f).astype("float32")
        observations.append(frame-np.median(frame))
    sky = np.median(observations, axis=0)
    scale = headers[0]["PIXSCALE"]
    anchor = np.array([158., 141.])  # zero-based x,y in first western-field image
    aligned = []
    for image, hdr in zip(observations, headers):
        residual = image - sky
        med = ndimage.median_filter(residual, 3)
        hot = np.abs(residual - med) > max(80, mad_std(residual)*7)
        residual[hot] = med[hot]
        filtered = ndimage.gaussian_filter(residual, .85)-ndimage.gaussian_filter(residual, 7)
        # At 270.7 deg rotator setting the detector +x,+y axes point south,west.
        dra = (hdr["RA"]-headers[0]["RA"]) * np.cos(np.deg2rad(headers[0]["DEC"]))
        ddec = hdr["DEC"]-headers[0]["DEC"]
        px, py = np.rint(anchor + np.array([ddec, dra])*3600/scale).astype(int)
        cut = filtered[py-12:py+13, px-12:px+13]
        dy, dx = np.array(np.unravel_index(cut.argmax(), cut.shape)) - 12
        measured = keck_centroid(filtered, (px+dx, py+dy))
        aligned.append(ndimage.shift(filtered, (anchor-measured)[::-1], order=1, cval=np.nan))
    return np.array(aligned), scale


def keck_rows(paths):
    aligned, scale = keck_stack(paths)
    rows = []
    for omit in [None, *range(len(aligned))]:
        sample = aligned if omit is None else np.delete(aligned, omit, axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            image = np.nanmedian(sample, axis=0)
        one = keck_centroid(image, (158,141))
        seven = keck_centroid(image, (103,228))
        east, north = -scale*(seven-one)[::-1]
        radius = 3.5
        ratio = (aperture(image, seven, radius, 1.8*radius, 2.8*radius) /
                 aperture(image, one, radius, 1.8*radius, 2.8*radius))
        rows.append([paths[omit].name if omit is not None else "none",east,north,ratio])
    return rows


def jw_centroid(image, xy, search_radius):
    x, y = np.rint(xy).astype(int)
    n = search_radius
    cut = image[y-n:y+n+1,x-n:x+n+1]
    # Bright control uses a smoothed peak. Use local high-pass for extended S2.
    smoothed = ndimage.gaussian_filter(np.nan_to_num(cut), 1 if n==9 else 1.4)
    if n != 9:
        smoothed -= ndimage.gaussian_filter(np.nan_to_num(cut),5)
    iy, ix = np.unravel_index(smoothed.argmax(), cut.shape)
    x, y = x + ix-n, y + iy-n
    patch = image[y-7:y+8,x-7:x+8]
    valid = np.isfinite(patch)
    dy, dx = np.mgrid[-7:8,-7:8]
    yy, xx, data = dy[valid], dx[valid], patch[valid]

    def residual(p):
        xc,yc,sx,sy,theta,amp,base,bx,by = p
        u = (xx-xc)*np.cos(theta) + (yy-yc)*np.sin(theta)
        v = -(xx-xc)*np.sin(theta) + (yy-yc)*np.cos(theta)
        return amp*np.exp(-.5*((u/sx)**2 + (v/sy)**2)) + base + bx*xx + by*yy - data

    background = np.median(data)
    initial = [0,0,1.5,1.3,0,max(data.max()-background,.001),background,0,0]
    bounds = ([-3,-3,.55,.55,-3.15,0,-np.inf,-np.inf,-np.inf],
              [3,3,6,6,3.15,np.inf,np.inf,np.inf,np.inf])
    fit = least_squares(residual, initial, bounds=bounds).x
    return np.array([x,y])+fit[:2], np.sort(fit[2:4])[::-1]*2.355


def jw_row(path):
    with fits.open(path) as hdus:
        image, header = hdus["SCI"].data, hdus["SCI"].header
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(header)
        filt = hdus[0].header["FILTER"]
        scale = np.sqrt(header["PIXAR_A2"])
        conversion = 1e12*header["PIXAR_SR"]  # MJy/sr to microJy per image pixel
        star, width = jw_centroid(image, wcs.all_world2pix(*REFERENCE_SKY,0),9)
        sky = np.array(wcs.all_pix2world(*star,0))
        radius = .22/scale
        fluxes = []
        for angle in np.linspace(0,2*np.pi,60,endpoint=False):
            rotation = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
            offset = rotation @ KECK_VECTOR
            p = wcs.all_world2pix(*(sky+offset/EAST_NORTH_SCALE),0)
            fluxes.append(conversion*aperture(image,p,radius,radius+1.5,1.65*radius+2))
        fluxes = np.array(fluxes)
        controls = fluxes[1:][np.isfinite(fluxes[1:])]
        position = ["","","",""]
        if filt == "F444W":
            seed = wcs.all_world2pix(*(sky+KECK_VECTOR/EAST_NORTH_SCALE),0)
            centroid, fwhm = jw_centroid(image,seed,12)
            east,north = (np.array(wcs.all_pix2world(*centroid,0))-sky)*EAST_NORTH_SCALE
            position = [east,north,fwhm[0]*scale,fwhm[1]*scale]
        return [filt,fluxes[0],np.median(controls),mad_std(controls),
                int(np.sum(controls>=fluxes[0])),len(controls),*position,width[0]*scale,width[1]*scale]


def write_csv(path, fields, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(fields.split())
        for row in rows:
            writer.writerow([f"{v:.5f}" if isinstance(v,(float,np.floating)) else v for v in row])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    with (HERE/"manifest.csv").open() as stream:
        manifest = list(csv.DictReader(stream))
    fetch(args.cache,manifest)
    args.output.mkdir(parents=True,exist_ok=True)
    keck = [args.cache/r["filename"] for r in manifest if r["instrument"]=="Keck/NIRC2"]
    jwst = [args.cache/r["filename"] for r in manifest if r["instrument"]=="JWST/NIRCam"]
    write_csv(args.output/"keck.csv","omitted_frame east_arcsec north_arcsec flux_ratio_7_to_1",keck_rows(keck))
    write_csv(args.output/"jwst.csv",
              "filter aperture_uJy rotated_median_uJy rotated_mad_sigma_uJy controls_ge_target controls_finite east_arcsec north_arcsec major_arcsec minor_arcsec anchor_major_arcsec anchor_minor_arcsec",
              [jw_row(f) for f in jwst])


if __name__ == "__main__":
    main()
