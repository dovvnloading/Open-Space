"""Verify the 2023 measurement in the separately processed VLASS continuum image."""

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from replay import Image, predicted

PRODUCT = (
    "VLASS3.1.se.T26t04.J045918+643000.06.2048.v1.I.iter3.image.pbcor.tt0.subim.fits"
)
URL = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/VLASS/" + PRODUCT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "image",
        type=Path,
        help="Original continuum image or a FITS cutout covering at least 40 arcsec around the target",
    )
    args = parser.parse_args()
    im = Image(args.image)
    pos = predicted(im.date)
    samples = np.array(
        [
            im.fit(pos, (x, y))
            for x in range(-50, 51, 10)
            for y in range(-50, 51, 10)
            if 15 < np.hypot(x, y) < 52
        ]
    )
    amp = im.fit(pos)
    med = np.median(samples)
    noise = np.median(abs(samples - med)) * 1.4826
    ra, dec, _, east, north = im.centroid(pos)
    print(
        json.dumps(
            {
                "parent_url": URL,
                "processed_image_sha256": hashlib.sha256(
                    args.image.read_bytes()
                ).hexdigest(),
                "date": im.date,
                "forced_mjy_beam": amp,
                "local_scatter_mjy_beam": noise,
                "local_control_count": len(samples),
                "local_control_median_mjy_beam": med,
                "centroid_ra": ra,
                "centroid_dec": dec,
                "offset_arcsec": np.hypot(east, north),
                "radio_beam_arcsec": im.axes.tolist(),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
