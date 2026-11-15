"""Shift and stack the 2023-08-14 CFHT images and the discovery-night controls."""
import argparse
import csv
import json
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, map_coordinates, shift

SERVICE = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync"
YY, XX = np.indices((141, 141)) - 70
PIXEL_SCALE = 0.1855  # tangent-plane grid, arcsec


def reproject(row, directory):
    product = row["product"]
    cache = directory / (product + ".fits")
    if not cache.exists():
        params = dict(ID=f"cadc:CFHT/{product}.fits.fz",
                      CIRCLE=f"{row['cut_ra_deg']} {row['cut_dec_deg']} 0.007")
        with urlopen(SERVICE + "?" + urlencode(params), timeout=120) as response:
            content = response.read()
        cache.write_bytes(content)
    ra, dec = float(row["fit_ra_deg"]), float(row["fit_dec_deg"])
    with fits.open(cache) as image:
        for chip in image:
            if chip.data is None:
                continue
            pixels = chip.data.astype(float)
            wcs = WCS(chip.header)
            x, y = wcs.all_world2pix([[ra, dec]], 0)[0]
            if not 8 < x < pixels.shape[1]-8 or not 8 < y < pixels.shape[0]-8:
                continue
            grid = np.column_stack((ra + XX.ravel()*PIXEL_SCALE/3600/np.cos(np.deg2rad(dec)),
                                    dec + YY.ravel()*PIXEL_SCALE/3600))
            px, py = wcs.all_world2pix(grid, 0).T
            stamp = map_coordinates(pixels, [py, px], order=1, mode="constant",
                                    cval=np.nan).reshape(XX.shape)
            median = np.nanmedian(stamp)
            return (stamp - median) / (1.4826 * np.nanmedian(abs(stamp - median)))
    return None  # CADC footprint includes a chip gap at this position.


def subtract_field(images, positions):
    offset = (positions - np.median(positions, axis=0)) * 3600/PIXEL_SCALE
    offset[:, 0] *= np.cos(np.deg2rad(positions[:, 1].mean()))
    result = []
    for i, image in enumerate(images):
        delta = offset - offset[i]
        selected = np.flatnonzero(np.linalg.norm(delta, axis=1) > 18)
        if len(selected) < 3:
            result.append(image)
            continue
        aligned = np.stack([shift(images[j], delta[j][::-1], order=1, mode="constant",
                                  cval=np.nan, prefilter=False) for j in selected])
        result.append(image - np.nanmedian(aligned, axis=0))
    return np.stack(result)


def significance(images):
    median = np.nanmedian(images, axis=0)
    mad = 1.4826 * np.nanmedian(abs(images-median), axis=0)
    clipped = np.where(abs(images-median) < np.maximum(3.3*mad, 1.), images, np.nan)
    combined = np.nanmean(clipped, axis=0)
    smooth = gaussian_filter(np.nan_to_num(combined), 2)
    ring = (XX**2 + YY**2 > 16**2) & (XX**2 + YY**2 < 52**2)
    base = np.median(smooth[ring])
    return (smooth - base) / (1.4826 * np.median(abs(smooth[ring] - base)))


def measure(rows, images):
    retained = [(row, image) for row, image in zip(rows, images) if image is not None]
    data = np.stack([image for _, image in retained])
    sky = np.array([[float(row["fit_ra_deg"]), float(row["fit_dec_deg"])] for row, _ in retained])
    residual = subtract_field(data, sky)
    sn = significance(residual)
    center = float(sn[70, 70])
    peak = float(np.nanmax(sn[XX**2+YY**2 <= 3**2]))
    comparisons = []
    for y in range(-42, 43, 9):
        for x in range(-42, 43, 9):
            if 15**2 < x*x+y*y < 47**2:
                comparisons.append(np.nanmax(sn[(XX-x)**2+(YY-y)**2 <= 3**2]))
    return dict(frames_available=len(rows), frames_on_detector=len(data),
                exposure_seconds=round(sum(float(row["exposure_s"]) for row, _ in retained), 2),
                center_snr=round(center, 2), raw_center_snr=round(significance(data)[70, 70], 2),
                half_center_snr=[round(significance(part)[70, 70], 2) for part in np.array_split(residual, 2)],
                peak_snr_in_0p56_arcsec=round(peak, 2), comparison_regions=len(comparisons),
                comparison_peaks_at_least_as_high=int(sum(value >= peak for value in comparisons)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, help="Keep cutouts in this directory between runs")
    args = parser.parse_args()
    rows = list(csv.DictReader((Path(__file__).resolve().parent / "inputs/precovery_frames.csv").open()))
    with tempfile.TemporaryDirectory() as scratch:
        directory = args.cache or Path(scratch)
        directory.mkdir(exist_ok=True, parents=True)
        for day in sorted(set(row["night_mjd"] for row in rows)):
            night = [row for row in rows if row["night_mjd"] == day]
            with ThreadPoolExecutor(max_workers=3) as pool:
                images = list(pool.map(lambda row: reproject(row, directory), night))
            print(json.dumps(dict(mjd=day, **measure(night, images))), flush=True)


if __name__ == "__main__":
    main()
