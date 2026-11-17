"""Forced local ePSF amplitudes on motion-propagated DR1 cutouts. Python 3, numpy, scipy, astropy."""
import argparse
import csv
import hashlib
import io
import math
import warnings
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy import log
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

log.setLevel("ERROR")


def stamp(table, header, x, y, xx, yy):
    """Evaluate the nearest detector ePSF at the native subpixel phase."""
    data = table.data
    xparent = x + 1 - header["CRPIX1A"]
    yparent = y + 1 - header["CRPIX2A"]
    k = np.argmin((data["XCENTER"] - xparent)**2 + (data["YCENTER"] - yparent)**2)
    grid = np.array(data["EPSF"][k], dtype=float)
    native_x = (grid.shape[1] - 1)/2 + (xx - x)*table.header["OVSMPX"]
    native_y = (grid.shape[0] - 1)/2 + (yy - y)*table.header["OVSMPY"]
    psf = map_coordinates(grid, [native_y, native_x], order=3, mode="constant", cval=0)
    return psf / np.sum(psf)


def fit_frame(row, payload, neighbors):
    with fits.open(io.BytesIO(payload), memmap=False) as data:
        im = data["IMAGE"].data.astype(float)
        var = data["VARIANCE"].data.astype(float)
        flags = data["FLAGS"].data
        h = data["IMAGE"].header
        wh = WCS(h)
        a, d = float(row["ra"]), float(row["dec"])
        x, y = map(float, wh.all_world2pix(a, d, 0))
        w = WCS(h, data, key="W")
        wave, bandwidth = map(float,w.wcs_pix2world(x, y, 0))
        if not (np.isfinite(x+y+wave) and 5 <= x < im.shape[1]-5 and 5 <= y < im.shape[0]-5):
            raise ValueError("edge_or_bad_wavelength")
        yy, xx = np.indices(im.shape)
        region = (xx-x)**2 + (yy-y)**2 <= 5**2
        allowed_flags = (1 << 21) | (1 << 12)  # SOURCE and FULLSAMPLE are informational.
        clean = (flags & ~allowed_flags) == 0
        valid = region & clean & np.isfinite(im) & np.isfinite(var) & (var > 0)
        central = (xx-x)**2 + (yy-y)**2 < 1.5**2
        if not np.all(clean[central]) or valid.sum() < 25:
            raise ValueError("quality_or_mask")
        target_psf = stamp(data["EPSF"], h, x, y, xx, yy)
        sky = [np.ones_like(xx), xx-x, yy-y]
        # Static catalog positions, explicit nuisance parameters, including close blends.
        refs = []
        for neighbor in neighbors:
            xn, yn = map(float, wh.all_world2pix(float(neighbor["ra"]),float(neighbor["dec"]),0))
            sep = math.hypot(xn-x, yn-y)
            if sep <= 5:
                refs.append((sep,stamp(data["EPSF"], h, xn, yn, xx, yy)))
        design = np.array([target_psf, *sky, *(psf for _,psf in refs)])[:,valid].T
        sigma = np.sqrt(var[valid])
        coef, _, rank, _ = np.linalg.lstsq(design/sigma[:,None], im[valid]/sigma, rcond=1e-9)
        if rank < len(coef):
            raise ValueError("degenerate_fit")
        inv = np.linalg.inv((design/sigma[:,None]).T @ (design/sigma[:,None]))
        chi2 = np.sum(((im[valid]-design@coef)/sigma)**2)/(valid.sum()-rank)
        # Local differential WCS, square degrees/pixel to sr/pixel, then MJy to microJy.
        rxy, dxy = wh.all_pix2world([x,x+0.1,x],[y,y,y+0.1],0)
        dra = (np.array(rxy[1:])-rxy[0]+180)%360-180
        jac = np.vstack([dra*math.cos(math.radians(dxy[0])), np.array(dxy[1:])-dxy[0]])/.1
        factor = abs(np.linalg.det(jac))*(math.pi/180)**2*1e12
        return dict(wavelength_um=wave,bandwidth_um=bandwidth,mjd=float(h["MJD-AVG"]),
                    flux_ujy=coef[0]*factor,sigma_ujy=math.sqrt(inv[0,0]*max(1,chi2))*factor,
                    chi2red=chi2,n_good=int(valid.sum()),n_neighbors=len(refs),
                    nearest_pixels=min((s for s,_ in refs),default=None))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--cache", type=Path)
    args = parser.parse_args()
    fields = ["target","cutout_url","ra","dec","sha256","status","wavelength_um","bandwidth_um",
              "mjd","flux_ujy","sigma_ujy","chi2red","n_good","n_neighbors","nearest_pixels"]
    with args.manifest.open() as inf, args.out.open("w") as outf:
        writer = csv.DictWriter(outf,fieldnames=fields)
        writer.writeheader()
        for row in csv.DictReader(inf):
            out = {k:row[k] for k in ("target","cutout_url","ra","dec")}
            try:
                local = args.cache/row["cache_path"] if args.cache else None
                if local and local.exists():
                    payload = local.read_bytes()
                else:
                    with urlopen(row["cutout_url"], timeout=90) as response:
                        payload = response.read()
                    if local:
                        local.parent.mkdir(parents=True,exist_ok=True)
                        local.write_bytes(payload)
                out["sha256"] = hashlib.sha256(payload).hexdigest()
                with (args.manifest.parent/row["neighbor_file"]).open() as file:
                    neighbors = list(csv.DictReader(file))
                out.update(fit_frame(row,payload,neighbors))
                out["status"] = "fit"
            except (OSError, ValueError, np.linalg.LinAlgError) as exc:
                out["status"] = str(exc)
            writer.writerow(out)


if __name__ == "__main__":
    main()
