"""Measure the fixed-position 73.8 MHz VLSSr brightness of ILTJ112022.08+502936.2."""

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

from astropy.io import fits
from photometry import RadioImage

SERVICE = "https://www.cv.nrao.edu/cgi-bin/newVLSSpostage.pl"
PARAMETERS = {
    "Equinox": "J2000",
    "RA": "11 20 22.08",
    "Dec": "+50 29 36.2",
    "Size": "1.0 1.0",
    "Cells": "15 15",
    "MAPROJ": "SIN",
    "rotate": "0",
    "Type": "application/octet-stream",
}


def measure(path):
    with tempfile.TemporaryDirectory() as tmp:
        adapted = Path(tmp) / "field.fits"
        with fits.open(path) as h:
            raw_sha = hashlib.sha256(h[0].data.tobytes()).hexdigest()
            hdr = h[0].header
            freq = hdr["CRVAL3"] / 1e6
            epoch = hdr["DATE-OBS"]
            if (
                abs(freq - 73.8) > 0.1
                or "VLSSR" not in " ".join(hdr["HISTORY"]).upper()
            ):
                raise ValueError("Unexpected VLSSr FITS product")
            # Postage service omits BMAJ/BMIN. VLSSr restored all images to 75".
            hdr["BMAJ"] = hdr["BMIN"] = 75 / 3600
            h.writeto(adapted)
        with RadioImage(adapted) as radio:
            result = radio.measure(170.092166, 50.493238)
    return dict(
        source="ILTJ112022.08+502936.2",
        field="1100+550",
        date_obs=epoch,
        freq_mhz=freq,
        raw_primary_pixels_sha256=raw_sha,
        clean_bias_correction_jy=0.66 * result["pixel_noise_jy_beam"],
        **result,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stamp", type=Path)
    p.add_argument("json_output", type=Path)
    a = p.parse_args()
    if not a.stamp.exists():
        cmd = ["curl", "-fsSL", "--retry", "2", "--max-time", "120"]
        for key, value in PARAMETERS.items():
            cmd.extend(["--data-urlencode", f"{key}={value}"])
        subprocess.run([*cmd, SERVICE, "-o", str(a.stamp)], check=True)
    result = measure(a.stamp)
    a.json_output.write_text(json.dumps(result, indent=2) + "\n")
