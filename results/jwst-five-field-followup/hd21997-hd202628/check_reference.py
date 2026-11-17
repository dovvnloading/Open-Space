"""Counterfactual RDI for both HD 202628 threshold crossings, same detector exposure."""
import argparse
import csv
from pathlib import Path
import sys

import numpy as np
from scipy import ndimage

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("merged", type=Path)
    p.add_argument("cache", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    sys.path.insert(0, str(a.merged / "nircam"))
    from nircam_reduce import cube, matched, median, offaxis_psf, reduce_cube, thresholds, write_table
    science, hdr, _ = cube(next(a.cache.glob("jw03989005001*long_calints.fits")))
    refs = {"HD190580": sorted(a.cache.glob("jw03989006001*long_calints.fits")),
            "HD21997_PSF": sorted(a.cache.glob("jw03989002001*long_calints.fits"))}
    xy = [hdr["CRPIX1"] - 1, hdr["CRPIX2"] - 1]
    scale = 3600*np.sqrt(abs(np.linalg.det([[hdr["CD1_1"],hdr["CD1_2"]],[hdr["CD2_1"],hdr["CD2_2"]]])))
    y,x = np.indices(science.shape[-2:])
    radius = np.hypot(x-xy[0], y-xy[1])*scale
    psf,_ = offaxis_psf(sorted(a.cache.glob("jw01537028*long_calints.fits"))[:2])
    kernel = psf[14:23,14:23]
    peaks = list(csv.DictReader((a.out / "HD202628_peaks.csv").open()))
    results = []
    for name, files in refs.items():
        data = [cube(f)[0] for f in files]
        residual = reduce_cube(science, data, radius)
        img = median(residual, axis=0)
        img[np.isfinite(residual).sum(0) < 4] = np.nan
        amp = matched(img, kernel)
        cut, snr,_ = thresholds(amp,radius,xy,scale)
        for n,s in enumerate(peaks,1):
            ix,iy = int(s["x0"]),int(s["y0"])
            sl = np.s_[iy-2:iy+3,ix-2:ix+3]
            rows = img[iy-4:iy+5,ix-4:ix+5]
            results.append(dict(id=n,x=ix,y=iy,reference=name,ref_star=cube(files[0])[2]["TARGNAME"],
                reference_dithers=len(files),pixel_MJy_sr=img[iy,ix],stamp_min_MJy_sr=np.nanmin(rows),
                stamp_max_MJy_sr=np.nanmax(rows),matched_snr=snr[iy,ix],max_snr_in_2px_box=np.nanmax(snr[sl]),
                n_crossings_in_2px_box=np.count_nonzero(amp[sl]>cut[sl])))
    write_table(a.out/"HD202628_reference_test.csv", list(results[0]), results)

if __name__=="__main__":
    main()
