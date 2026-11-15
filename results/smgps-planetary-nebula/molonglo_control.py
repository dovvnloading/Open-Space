"""Compare 843 MHz and 1.3 GHz on the native Molonglo grid at a common beam."""

import csv
import hashlib
import json
import math
from urllib.request import urlopen

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from radio_photometry import HERE, POSITION, RADIO, RADIO_SHA
from scipy.ndimage import gaussian_filter, map_coordinates

URL = "https://skyview.gsfc.nasa.gov/surveys/sumss/mosaics/Galactic/J0936M52.FITS"
SHA = "e30f277927291064195e232fdbd75effb8ac33f6936158d4a5661c82a6537df5"


def main():
    path = HERE / "cache" / "J0936M52.FITS"
    path.parent.mkdir(exist_ok=True)
    if not path.exists():
        with urlopen(URL, timeout=90) as response:
            path.write_bytes(response.read())
    for name, checksum in [(path, SHA), (RADIO, RADIO_SHA)]:
        if hashlib.sha256(name.read_bytes()).hexdigest() != checksum:
            raise ValueError(f"Unexpected input: {name}")
    with fits.open(path) as h:
        mol, hdr = h[0].data.copy(), h[0].header
    with fits.open(RADIO) as h:
        meer, mw = h[0].data.copy(), WCS(h[0].header, naxis=2)
    w = WCS(hdr, naxis=2)
    if hdr["BPA"] != 0 or hdr["BUNIT"].strip().upper() != "JY/BEAM":
        raise ValueError("Unexpected Molonglo beam or units")
    maj, minr = hdr["BMAJ"] * 3600, hdr["BMIN"] * 3600
    cx, cy = w.world_to_pixel(POSITION)
    ys, xs = np.mgrid[
        round(float(cy)) - 40 : round(float(cy)) + 41,
        round(float(cx)) - 40 : round(float(cx)) + 41,
    ]
    sky = w.pixel_to_world(xs, ys).icrs
    sky = SkyCoord(ra=sky.ra, dec=sky.dec, frame="icrs")
    east, north = (v.arcsec for v in POSITION.spherical_offsets_to(sky))
    r = np.hypot(east, north)
    az = (np.degrees(np.arctan2(east, north)) + 360) % 360
    # Reproject the 8" data to a fine equatorial grid before convolution.
    # Keep the zero-padded boundary stress separate below; trust the central aperture only.
    grid = np.arange(-210, 211, 1.5)
    eg, ng = np.meshgrid(grid, grid)
    coords = POSITION.spherical_offsets_by(eg * u.arcsec, ng * u.arcsec)
    mx, my = mw.world_to_pixel(coords)
    arr = map_coordinates(meer, [my, mx], order=1, mode="constant", cval=0)
    sigma = np.sqrt(np.array([maj, minr]) ** 2 - 8**2) / (
        2 * math.sqrt(2 * math.log(2)) * 1.5
    )
    smooth = gaussian_filter(arr, sigma) * maj * minr / 8**2
    # Same result sampled onto the 843 MHz grid. Jy per common restoring beam.
    smeared = map_coordinates(
        smooth,
        [(north - grid[0]) / 1.5, (east - grid[0]) / 1.5],
        order=1,
        mode="constant",
        cval=np.nan,
    )
    pixel_area = abs(hdr["CDELT1"] * hdr["CDELT2"]) * 3600**2
    k = 1000 * pixel_area / (math.pi / (4 * math.log(2)) * maj * minr)
    # Compare normalized convolution at the finite MeerKAT cutout boundary.
    valid = (
        (mx >= 0) & (mx <= meer.shape[1] - 1) & (my >= 0) & (my <= meer.shape[0] - 1)
    )
    weight = gaussian_filter(valid.astype(float), sigma)
    corrected = np.divide(
        smooth, weight, out=np.zeros_like(smooth), where=weight > 0.05
    )
    alt = map_coordinates(
        corrected,
        [(north - grid[0]) / 1.5, (east - grid[0]) / 1.5],
        order=1,
        mode="constant",
        cval=np.nan,
    )
    images = {"MGPS_843_MHz": mol[ys, xs], "SMGPS_1300_MHz_matched": smeared}
    rows = []
    for label, im in images.items():
        for ap in [85, 90, 100, 110]:
            mask = r <= ap
            for ri, ro in [(110, 135), (120, 150), (135, 165), (150, 175)]:
                ann = (r > ri) & (r < ro)
                # Median annuli reduce contamination by the bright neighbour.
                med = float(np.nanmedian(im[ann]))
                rows.append(
                    [
                        label,
                        ap,
                        ri,
                        ro,
                        "all",
                        float(np.sum(im[mask] - med) * k),
                        med * 1000,
                    ]
                )
        for j in range(8):
            ann = (r > 120) & (r < 150) & (az >= j * 45) & (az < (j + 1) * 45)
            med = float(np.nanmedian(im[ann]))
            rows.append(
                [
                    label,
                    90,
                    120,
                    150,
                    str(j),
                    float(np.sum(im[r <= 90] - med) * k),
                    med * 1000,
                ]
            )
    with (HERE / "molonglo-apertures.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            [
                "image",
                "radius_arcsec",
                "sky_inner_arcsec",
                "sky_outer_arcsec",
                "octant_east_of_north",
                "flux_mjy",
                "sky_mjy_per_beam",
            ]
        )
        writer.writerows(rows)
    summary = {}
    for label, im in images.items():
        subset = [q for q in rows if q[0] == label and q[1] == 90]
        full = [q[5] for q in subset if q[4] == "all"]
        sector = [q[5] for q in subset if q[4] != "all"]
        neighbour = np.where((r > 120) & (r < 240), im, -np.inf)
        # Molonglo source outside the nebular aperture; its foreground/background
        # status is unknown. Do not subtract it from the nebula.
        j = np.unravel_index(np.nanargmax(neighbour), im.shape)
        summary[label] = {
            "standard_mjy": full[1],
            "annulus_range_mjy": [min(full), max(full)],
            "sky_octant_range_mjy": [min(sector), max(sector)],
            "brightest_nearby_peak_mjy_per_beam": float(im[j] * 1000),
            "nearby_peak_offset_arcsec": float(r[j]),
        }
    # Nearby blank-sky controls include the Galactic-plane sky and any unrelated emitters.
    # Use the same on-source 90" aperture, annulus 120-150", and no clipping.
    blank = []
    for angle in range(0, 360, 45):
        rad = math.radians(angle)
        rr = np.hypot(east - 285 * math.sin(rad), north - 285 * math.cos(rad))
        ann = (rr > 120) & (rr < 150)
        bg = float(np.median(images["MGPS_843_MHz"][ann]))
        blank.append(float(np.sum(images["MGPS_843_MHz"][rr <= 90] - bg) * k))
    summary.update(
        {
            "matched_smoothed_edge_normalization_check_mjy": float(
                np.sum(alt[r <= 90] - np.nanmedian(alt[(r > 120) & (r < 150)])) * k
            ),
            "mosaic_sha256": SHA,
            "mosaic_url": URL,
            "common_beam_major_minor_arcsec": [maj, minr],
            "approx_molonglo_pixel_area_sq_arcsec": pixel_area,
            "mgps_eight_offsource_aperture_mjy_at_radius_285_arcsec": blank,
        }
    )
    (HERE / "molonglo-control.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
