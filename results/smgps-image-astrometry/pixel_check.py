"""Validate native images against an independent radio survey, with shifted controls."""

import csv
import hashlib
from urllib.parse import urlencode
from urllib.request import urlopen

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from harvest import HERE, filename_url
from scipy.ndimage import maximum_filter
from scipy.spatial import cKDTree

FILES = [
    "SMGPS_Candidate_G040.449+00.540.fits",
    "SMGPS_Candidate_G043.942+01.442.fits",
    "SMGPS_Candidate_G046.577+00.221.fits",
    "SMGPS_Candidate_G050.489-00.410.fits",
    "SMGPS_Candidate_G051.746+01.384.fits",
    "SMGPS_Candidate_G057.868+00.221.fits",
]
NVSS = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"


def obtain(url, dest):
    if not dest.exists():
        with urlopen(url, timeout=90) as resp:
            data = resp.read()
        dest.write_bytes(data)
    return dest


def main():
    cache = HERE / "pixel-cache"
    cache.mkdir(exist_ok=True)
    controls = [(0, 0), (50, 0), (-50, 0), (0, 50), (0, -50)]
    matches = []
    with (HERE / "pixel-check.csv").open("w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(
            [
                "filename",
                "image_sha256",
                "nvss_sha256",
                "markers",
                "native_within_6as",
                "control_px_plus50_x",
                "control_px_minus50_x",
                "control_px_plus50_y",
                "control_px_minus50_y",
                "native_within_12as",
                "control_12as_plus50_x",
                "control_12as_minus50_x",
                "control_12as_plus50_y",
                "control_12as_minus50_y",
                "median_match_arcsec",
            ]
        )
        for file in FILES:
            image_path = obtain(filename_url(file), cache / file)
            with fits.open(image_path) as hdus:
                wcs = WCS(hdus[0].header, naxis=2)
                im = hdus[0].data
                ny, nx = im.shape
                center = wcs.pixel_to_world((nx - 1) / 2, (ny - 1) / 2).icrs
            args = {
                "-source": "VIII/65/nvss",
                "-c": f"{center.ra.deg:.6f} {center.dec.deg:+.6f}",
                "-c.rm": f"{0.0125 * min(nx, ny) + 6:.5f}",
                "-out": "NVSS,RAJ2000,DEJ2000,S1.4,MajAxis,MinAxis",
                "-out.max": "unlimited",
                "-oc.form": "dec",
            }
            nvss_path = obtain(
                NVSS + urlencode(args), (cache / file).with_suffix(".tsv")
            )
            # Hash the data records, excluding the request timestamp in the preamble.
            lines = [
                s
                for s in nvss_path.read_text().splitlines()
                if len(s.split("\t")) == 6 and s[:6].isdigit()
            ]
            if "VizieR Astronomical Server" not in nvss_path.read_text() or not lines:
                raise ValueError(f"Unverified NVSS response for {file}")
            ra, dec, ids = [], [], []
            for line in lines:
                name, alpha, delta, flux, major, _minor = line.strip().split("\t")
                if float(flux) >= 10 and float(major.lstrip("< ")) <= 45:
                    ids.append(name.strip())
                    ra.append(alpha.strip())
                    dec.append(delta.strip())
            ref = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame="icrs").galactic
            xy = np.column_stack(wcs.world_to_pixel(ref))
            keep = np.all((xy >= 55) & (xy <= np.array([nx - 1, ny - 1]) - 55), axis=1)
            xy, ids = xy[keep], np.asarray(ids)[keep]
            # 1.5" per pixel. Peaks must be local maxima in a 16.5" box and pass 8 robust sigma.
            med = np.nanmedian(im)
            rms = 1.4826 * np.nanmedian(abs(im - med))
            clean = np.nan_to_num(im, nan=-np.inf)
            peaks = np.argwhere(
                (clean > med + 8 * rms) & (clean == maximum_filter(clean, size=11))
            )[:, ::-1]
            tree = cKDTree(peaks)
            ds = [tree.query(xy + delta)[0] * 1.5 for delta in controls]
            for name, pos, d in zip(ids, xy, ds[0]):
                matches.append([file, name, *[f"{x:.4f}" for x in pos], f"{d:.4f}"])
            counts = [int(np.sum(ds_i <= 6)) for ds_i in ds]
            med_match = np.median(ds[0][ds[0] <= 6]) if counts[0] else np.nan
            writer.writerow(
                [
                    file,
                    hashlib.sha256(image_path.read_bytes()).hexdigest(),
                    hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest(),
                    len(ids),
                    *counts,
                    *[int(np.sum(d <= 12)) for d in ds],
                    f"{med_match:.3f}",
                ]
            )
    with (HERE / "nvss-markers.csv").open("w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(
            [
                "filename",
                "nvss",
                "predicted_x_zero_based",
                "predicted_y_zero_based",
                "nearest_radio_peak_arcsec",
            ]
        )
        writer.writerows(matches)


if __name__ == "__main__":
    main()
