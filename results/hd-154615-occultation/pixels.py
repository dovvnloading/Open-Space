"""Rebuild TESS eclipse difference images and the secondary sightline's harmonic map."""
import argparse
import csv
import hashlib
import io
from pathlib import Path
import urllib.parse
import urllib.request
import zipfile

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
HD = "5965875226751271040"
ATO = "4162228724694601600"
SECTORS = ((HD, 39), (HD, 93), (ATO, 92))


def difference(t, flux, phase):
    folded = (t - 1644.95) / 1.851076 - phase
    cycle = np.floor(folded + 0.5)
    distance = abs(folded - cycle)
    events = []
    for number in np.unique(cycle):
        same = cycle == number
        inside = same & (distance < 0.026)
        outside = same & (distance > 0.07) & (distance < 0.16)
        if inside.sum() >= 2 and outside.sum() >= 5:
            events.append(np.nanmedian(flux[outside], axis=0) - np.nanmedian(flux[inside], axis=0))
    return np.nanmedian(events, axis=0), len(events)


def harmonic_map(t, flux, xy):
    phase = 2 * np.pi * t / (1.43860334 / 2)
    block = np.floor((t - t.min()) / 1.5)
    estimates = []
    for number in np.unique(block):
        sel = block == number
        local_time = t[sel] - np.mean(t[sel])
        if np.ptp(local_time) < 0.8:
            continue
        design = np.column_stack((np.ones(sel.sum()), local_time, local_time**2,
                                  np.sin(phase[sel]), np.cos(phase[sel])))
        coefs = np.linalg.lstsq(design, flux[sel].reshape(sel.sum(), 81), rcond=None)[0]
        estimates.append(coefs[-2:].reshape(2, 9, 9))
    sine_cosine = np.median(estimates, axis=0)
    x, y = np.rint(xy).astype(int)
    direction = sine_cosine[:, y, x]
    return np.sum(direction[:, None, None] * sine_cosine, axis=0) / np.linalg.norm(direction), len(estimates)


def fit_position(image):
    y, x = np.indices(image.shape)
    image = image / np.max(image)
    def residual(params):
        cx, cy, sx, sy, theta, height, bg = params
        rx = (x-cx)*np.cos(theta) + (y-cy)*np.sin(theta)
        ry = (y-cy)*np.cos(theta) - (x-cx)*np.sin(theta)
        return (bg + height*np.exp(-0.5*((rx/sx)**2+(ry/sy)**2)) - image).ravel()
    py, px = np.unravel_index(np.argmax(image), image.shape)
    fit = least_squares(residual, [px, py, 0.9, 0.9, 0, 1, 0],
                        bounds=([-1, -1, 0.35, 0.35, -4, 0, -1], [9, 9, 3, 3, 4, 2, 1]),
                        loss="soft_l1", f_scale=0.06)
    return fit.x[:2]


def main(cache):
    cache.mkdir(parents=True, exist_ok=True)
    with (ROOT / "inputs" / "gaia_sources.csv").open() as stream:
        stars = {row["source_id"]: row for row in csv.DictReader(stream)}
    summary, maps = [], []
    for source, sector in SECTORS:
        star = stars[source]
        args = dict(ra=star["ra"], dec=star["dec"], x=9, y=9, units="px", sector=sector)
        url = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?" + urllib.parse.urlencode(args)
        path = cache / f"pixels_{source}_{sector}.zip"
        if not path.exists():
            path.write_bytes(urllib.request.urlopen(url, timeout=180).read())
        with zipfile.ZipFile(path) as archive:
            contents = archive.read(next(n for n in archive.namelist() if n.endswith(".fits")))
        with fits.open(io.BytesIO(contents)) as hdus:
            data = hdus[1].data
            xy = WCS(hdus[2].header).all_world2pix([[float(star["ra"]), float(star["dec"])]], 0)[0]
            good = (data["QUALITY"] == 0) & np.isfinite(data["TIME"])
            t, flux = data["TIME"][good], data["FLUX"][good]
            pixel_hash = hashlib.sha256(np.asarray(flux, dtype=">f4").tobytes()).hexdigest()
            if source == HD:
                selections = [(name, difference(t, flux, phase)) for name, phase in
                              (("primary", 0), ("secondary", 0.4975))]
            else:
                selections = [("twice_orbital_frequency", harmonic_map(t, flux, xy))]
            for label, (image, count) in selections:
                cx, cy = fit_position(image)
                summary.append(dict(source_id=source, sector=sector, component=label, events_or_blocks=count,
                                    target_x=xy[0], target_y=xy[1], fitted_x=cx, fitted_y=cy,
                                    offset_pixels=np.hypot(cx-xy[0], cy-xy[1]), cadence_count=len(t),
                                    flux_sha256=pixel_hash, url=url))
                for (y, x), amplitude in np.ndenumerate(image):
                    maps.append(dict(source_id=source, sector=sector, component=label, x=x, y=y,
                                     amplitude_electrons_per_s=amplitude))
    for file, data in [("localization.csv", summary), ("difference_images.csv", maps)]:
        with (ROOT / file).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, list(data[0]))
            writer.writeheader()
            writer.writerows(data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    main(parser.parse_args().cache)
