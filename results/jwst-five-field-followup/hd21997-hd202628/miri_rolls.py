"""Split-roll, detector-fixed and alternate-reference tests for the MIRI peaks."""

import argparse
import csv
from pathlib import Path
import sys
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import minimum_filter


def read(f):
    with fits.open(f) as hd:
        assert hd[0].header["FILTER"] == "F1140C"
        z = hd["SCI"].data.astype(float)
        z[(hd["DQ"].data & 3) != 0] = np.nan
        return z, hd["SCI"].header.copy(), hd[0].header.copy()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("merged", type=Path)
    p.add_argument("cache", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    sys.path[:0] = [str(a.merged / "nircam"), str(a.merged / "miri")]
    from nircam_reduce import reduce_cube, median, write_table
    from search import smooth_response

    rows = list(csv.DictReader((a.out / "auxiliary_peaks.csv").open()))
    res = []
    for path in sorted(a.cache.glob("*mirimage_c10*_psfsub.fits")):
        pipeline, h, phead = read(path)
        field = phead["TARGNAME"].replace(" ", "")
        w = WCS(h, naxis=2)
        cx, cy = h["CRPIX1"] - 1, h["CRPIX2"] - 1
        center = w.pixel_to_world(cx, cy)
        yy, xx = np.indices(pipeline.shape[-2:])
        scale = np.sqrt(h["PIXAR_SR"]) * 206264.806
        rad = np.hypot(xx - cx, yy - cy) * scale
        native, _, _ = read(
            a.cache
            / path.name.replace("_" + path.name.split("_")[-2] + "_psfsub", "_calints")
        )
        stacks = [("pipeline", pipeline)]
        for prefix in ["jw02538009", "jw02538014"]:
            files = sorted(a.cache.glob(prefix + "*mirimage_calints.fits"))
            assert len(files) == 9
            refs = [read(f)[0] for f in files]
            star = read(files[0])[2]["TARGNAME"]
            stacks.append(("RDI_" + star, reduce_cube(native, refs, rad)))
        for method, stack in stacks:
            for subset, arr in [
                ("all", stack),
                ("odd", stack[::2]),
                ("even", stack[1::2]),
            ]:
                im = median(arr, axis=0)
                valid = np.isfinite(im)
                good = (
                    minimum_filter(valid, size=5, mode="constant", cval=0)
                    & (abs(xx - cx) * scale > 0.37)
                    & (abs(yy - cy) * scale > 0.37)
                )
                filt = smooth_response(np.where(valid, im, median(im)))
                for s in [
                    r for r in rows if r["field"] == field and r["band"] == "F1140C"
                ]:
                    for loc in ["sky", "detector"]:
                        x, y = float(s["x"]), float(s["y"])
                        if loc == "sky":
                            en = np.array([s["east"], s["north"]], float)
                            x, y = w.world_to_pixel(
                                center.spherical_offsets_by(*(en * u.arcsec))
                            )
                            x, y = float(x), float(y)
                        ix, iy = int(round(x)), int(round(y))
                        r0 = rad[iy, ix]
                        lo = 2.5 + np.floor((r0 - 2.5) / 0.5) * 0.5
                        ring = good & (rad >= lo) & (rad < lo + 0.5)
                        m = median(filt[ring])
                        sig = 1.4826 * median(abs(filt[ring] - m))
                        res.append(
                            dict(
                                field=field,
                                original=s["file"],
                                id=s["id"],
                                roll=path.name,
                                method=method,
                                subset=subset,
                                track=loc,
                                x=x,
                                y=y,
                                usable=int(good[iy, ix]),
                                fixed_snr=(filt[iy, ix] - m) / sig,
                            )
                        )
    write_table(a.out / "MIRI_roll_tests.csv", list(res[0]), res)


if __name__ == "__main__":
    main()
