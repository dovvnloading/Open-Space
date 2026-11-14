"""Fit the planetary transit and resolved binary variability to every TESS pixel."""

import csv

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def locate(path, field_path, binary_flux, planet_epoch, planet_period, transit_shape):
    with field_path.open() as handle:
        field = {row["source_id"]: row for row in csv.DictReader(handle)}
    host = field["4534144923690481408"]
    neighbor = field["4534144923690486144"]
    with fits.open(path) as product:
        data = product[1].data
        cube = data["FLUX"].astype(float)
        good = ((data["QUALITY"] == 0) & np.isfinite(data["TIME"])
                & np.all(np.isfinite(cube), axis=(1, 2)))
        time = data["TIME"][good]
        spectra = cube[good].reshape((good.sum(), -1))
        sector = product[0].header["SECTOR"]
        wcs = WCS(product[2].header)
        height, width = cube.shape[1:]
    offset = (time - planet_epoch + planet_period / 2) % planet_period - planet_period / 2
    signals = np.column_stack((binary_flux(time), -transit_shape(offset)))
    full = np.column_stack((signals, spectra))
    day = np.floor(time)
    chunks = []
    for date in np.unique(day):
        use = day == date
        if use.sum() < 100:
            continue
        design = np.column_stack((np.ones(use.sum()), time[use] - date))
        coef = np.linalg.lstsq(design, full[use], rcond=None)[0]
        chunks.append(full[use] - design @ coef)
    residual = np.vstack(chunks)
    maps = np.linalg.lstsq(residual[:, :2], residual[:, 2:], rcond=None)[0]
    yy, xx = np.indices((height, width))
    output = []
    for name, frame in zip(("neighbor_waveform", "planet_transit"), maps.reshape(2, height, width)):
        peak_y, peak_x = np.unravel_index(frame.argmax(), frame.shape)
        for radius in (1, 2, 3):
            cutout = (abs(yy - peak_y) <= radius) & (abs(xx - peak_x) <= radius)
            background = np.median(frame[~cutout])
            weights = np.clip(frame - background, 0, None) * cutout
            x, y = (np.sum(weights * coord) / weights.sum() for coord in (xx, yy))
            ra, dec = wcs.pixel_to_world_values(x, y)
            # Gaia J2016.0 positions; propagated motion to 2024 is under 0.3 arcsec.
            def distance_arcsec(source):
                first = (float(source["ra"]) - ra) * np.cos(np.deg2rad(dec))
                return 3600 * np.hypot(first, float(source["dec"]) - dec)

            output.append(dict(sector=sector, signal=name, width_pixels=2 * radius + 1,
                               centroid_x=x, centroid_y=y, ra_deg=ra, dec_deg=dec,
                               host_arcsec=distance_arcsec(host),
                               neighbor_arcsec=distance_arcsec(neighbor)))
    return output
