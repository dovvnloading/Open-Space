"""Tie the later ALMA image to the published S2–S4 centroids and test S1."""

import csv
import io
import sys
import tarfile
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from scipy.optimize import least_squares
from auxiliary import apparent


def main():
    cache, out = map(Path, sys.argv[1:3])
    here = Path(__file__).resolve().parent
    catalog = list(csv.DictReader((here / "catalogs" / "HD202628_SIMBAD.csv").open()))
    host = next(
        r
        for r in csv.DictReader((here / "catalogs" / "HD202628_GaiaDR3.csv").open())
        if r["source_id"] == "6577318075351824512"
    )
    with tarfile.open(cache / "alma2018.tar") as t:
        n = next(n for n in t.getnames() if "HD202628_sci" in n and n.endswith("fits"))
        hd = fits.open(io.BytesIO(t.extractfile(n).read()))
    h = hd[0].header
    data = hd[0].data.squeeze().astype(float) * 1e6
    w = WCS(h).celestial
    scale = w.proj_plane_pixel_scales()[0].value * 3600
    dd = apparent(host, Time(h["DATE-OBS"]).mjd) - apparent(
        host, Time("2017-04-22").mjd
    )
    result = []
    for i in [1, 2, 3, 4]:
        r = next(r for r in catalog if r["main_id"].endswith("S" + str(i)))
        sky = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
        x0, y0 = w.world_to_pixel(sky)
        yy, xx = np.indices(data.shape)
        dx, dy = (xx - x0) * scale, (yy - y0) * scale
        for radius in [1.1, 1.4, 1.7]:
            mask = (np.hypot(dx, dy) < radius) & np.isfinite(data)
            x = dx[mask]
            y = dy[mask]
            z = data[mask]

            def residual(p):
                amp, cx, cy, sx, sy, rho, b0, bx, by = p
                u, v = (x - cx) / sx, (y - cy) / sy
                return (
                    b0
                    + bx * x
                    + by * y
                    + amp
                    * np.exp(-0.5 * (u * u + v * v - 2 * rho * u * v) / (1 - rho * rho))
                    - z
                )

            sol = least_squares(
                residual,
                [np.max(z) - np.median(z), 0, 0, 0.45, 0.45, 0, np.median(z), 0, 0],
                bounds=(
                    [0, -0.7, -0.7, 0.2, 0.2, -0.8, -np.inf, -np.inf, -np.inf],
                    [np.inf, 0.7, 0.7, 1.5, 1.5, 0.8, np.inf, np.inf, np.inf],
                ),
            )
            cen = w.pixel_to_world(x0 + sol.x[1] / scale, y0 + sol.x[2] / scale)
            east, north = sky.spherical_offsets_to(cen)
            result.append(
                dict(
                    source="S" + str(i),
                    fit_radius_arcsec=radius,
                    shift_east=east.arcsec,
                    shift_north=north.arcsec,
                )
            )
    for r in result:
        refs = [
            z
            for z in result
            if z["fit_radius_arcsec"] == r["fit_radius_arcsec"] and z["source"] != "S1"
        ]
        offsets = np.array([[z["shift_east"], z["shift_north"]] for z in refs])
        tie = np.median(offsets, axis=0)
        de = np.array([r["shift_east"], r["shift_north"]]) - tie
        r.update(
            tie_east=tie[0],
            tie_north=tie[1],
            tie_rms=np.sqrt(np.mean(np.sum((offsets - tie) ** 2, axis=1))),
            corrected_east=de[0],
            corrected_north=de[1],
            star_motion_east=dd[0],
            star_motion_north=dd[1],
            resid_stationary=np.linalg.norm(de),
            resid_comoving=np.linalg.norm(de - dd),
        )
    with (out / "ALMA_S1_motion.csv").open("w") as fd:
        writer = csv.DictWriter(fd, fieldnames=list(result[0]))
        writer.writeheader()
        writer.writerows(result)


if __name__ == "__main__":
    main()
