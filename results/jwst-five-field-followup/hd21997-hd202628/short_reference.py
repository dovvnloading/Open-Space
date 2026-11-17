"""Test every F200W mosaic peak in detector integrations and two disjoint RDI libraries."""

import argparse
import csv
from pathlib import Path
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("merged", type=Path)
    p.add_argument("cache", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    sys.path[:0] = [str(a.merged / "nircam"), str(a.merged / "miri")]
    from nircam_check_sources import geometry, read_short
    from nircam_reduce import median, reduce_cube, write_table
    from search import smooth_response

    peaks = list(csv.DictReader((a.out / "auxiliary_peaks.csv").open()))
    results = []
    for field, obs in [("HD21997", "001"), ("HD202628", "005")]:
        science_path = next(a.cache.glob("jw03989" + obs + "*nrca2_calints.fits"))
        sci, h = read_short(science_path)
        w = WCS(h, naxis=2)
        cd, cen, scale, r = geometry(h, sci.shape[1:])
        primary = fits.getheader(science_path)
        host = w.pixel_to_world(
            *w.world_to_pixel_values(primary["TARG_RA"], primary["TARG_DEC"])
        )
        libs = [("HD25945", "002"), ("HD190580", "006")]
        products = [("science", sci)]
        for star, n in libs:
            refs = [
                read_short(f)[0]
                for f in sorted(a.cache.glob("jw03989" + n + "*nrca2_calints.fits"))
            ]
            products.append((star, reduce_cube(sci, refs, r, sky_bounds=(3.4, 4.2))))
        for library, stack in products:
            image = median(stack, axis=0)
            good = np.isfinite(image)
            filt = smooth_response(np.where(good, image, np.nanmedian(image)))
            for row in peaks:
                if row["field"] != field or row["band"] != "F200W":
                    continue
                en = np.array([row["east"], row["north"]], float)
                src = host.spherical_offsets_by(*(en * u.arcsec))
                x, y = w.world_to_pixel(src)
                ix, iy = int(round(float(x))), int(round(float(y)))
                ring = good & (np.abs(r - np.hypot(*en)) < 0.25)
                scatter = 1.4826 * median(abs(filt[ring] - median(filt[ring])))
                patch = filt[iy - 4 : iy + 5, ix - 4 : ix + 5]
                peak = np.unravel_index(np.argmax(patch), patch.shape)
                score = (np.max(patch) - median(filt[ring])) / scatter
                results.append(
                    dict(
                        field=field,
                        id=row["id"],
                        reference=library,
                        x_predicted=x,
                        y_predicted=y,
                        max_in_4pix_snr=score,
                        peak_x=ix + peak[1] - 4,
                        peak_y=iy + peak[0] - 4,
                        filtered_value_at_predicted=filt[iy, ix],
                        native_core_median=median(
                            image[iy - 1 : iy + 2, ix - 1 : ix + 2]
                        ),
                        pixel_MJy_sr=image[iy, ix],
                        ring_sigma=scatter,
                    )
                )
    write_table(a.out / "F200W_reference_test.csv", list(results[0]), results)


if __name__ == "__main__":
    main()
