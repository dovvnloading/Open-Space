"""Download independent July 2024 CADC cutouts and measure local aperture depth."""

import hashlib
import io
import json
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve

SERVICE = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync"
# Cutout center is fixed to preserve the requested product. Test coordinates are
# at the midpoint of the exposure. They can also be recomputed with predict.py.
FRAMES = [
    ("3030533p", 350.68068, -6.50420, 350.68063733, -6.50420954),
    ("3030535p", 349.17030, -6.78115, 349.17024302, -6.78117519),
]


def inspect(image, ra, dec):
    for chip in image:
        if chip.data is None:
            continue
        pixels = chip.data.astype(float)
        ysize, xsize = pixels.shape
        x, y = WCS(chip.header).all_world2pix([[ra, dec]], 0)[0]
        if not (16 < x < xsize - 16 and 16 < y < ysize - 16):
            continue
        yy, xx = np.indices((33, 33)) - 16
        radius2 = xx**2 + yy**2
        aperture = radius2 <= 25
        annulus = (radius2 >= 100) & (radius2 <= 256)
        kernel = aperture.astype(float) - annulus * (aperture.sum() / annulus.sum())
        pixels[~np.isfinite(pixels)] = np.nanmedian(pixels)
        flux = fftconvolve(pixels, kernel, mode="same")
        gy, gx = np.indices(pixels.shape)
        nearby = (np.hypot(gx - x, gy - y) < 240) & (gx > 18) & (gy > 18)
        nearby &= (gx < xsize - 18) & (gy < ysize - 18)
        sample = flux[::12, ::12][nearby[::12, ::12]]
        background = np.median(sample)
        noise = 1.4826 * np.median(np.abs(sample - background))
        significance = (flux[round(y), round(x)] - background) / noise
        hdr = chip.header
        zp = hdr["PHOT_C"] + hdr["PHOT_K"] * (hdr["AIRMASS"] - 1)
        return dict(chip=chip.name, x_zero_based=x, y_zero_based=y,
                    samples=len(sample), sigma_adu=noise, target_snr=significance,
                    uncorr_aperture_mag_5sigma=zp - 2.5*np.log10(5 * noise / hdr["EXPTIME"]))
    raise ValueError("The expected direction fell outside every covered detector")


def main():
    for product, cx, cy, ra, dec in FRAMES:
        query = urlencode(dict(ID=f"cadc:CFHT/{product}.fits.fz", CIRCLE=f"{cx} {cy} 0.035"))
        with urlopen(SERVICE + "?" + query, timeout=120) as response:
            content = response.read()
        with fits.open(io.BytesIO(content)) as image:
            info = inspect(image, ra, dec)
        print(json.dumps(dict(product=product, fits_sha256=hashlib.sha256(content).hexdigest(), **info)))


if __name__ == "__main__":
    main()
