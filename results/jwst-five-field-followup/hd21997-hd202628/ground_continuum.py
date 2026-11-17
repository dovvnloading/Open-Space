"""Test stationary and host-following positions in archived ALMA continuum maps."""

import argparse
import csv
import io
import sys
import tarfile
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
from scipy.ndimage import map_coordinates

from auxiliary import GAIA, apparent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("merged", type=Path)
    p.add_argument("cache", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    sys.path.insert(0, str(a.merged / "nircam"))
    from nircam_reduce import write_table

    here = Path(__file__).resolve().parent
    cats = list(csv.DictReader((a.out / "catalog_matches.csv").open()))
    peaks = list(csv.DictReader((a.out / "auxiliary_peaks.csv").open()))
    results = []
    for field, pkg, pattern in [
        (
            "HD202628",
            "alma2018.tar",
            "HD202628_sci.spw19_21_23_25.mfs.I.manual.image.pbcor.fits",
        ),
        ("HD21997", "alma21997.tar", "HD21997_sci"),
    ]:
        with tarfile.open(a.cache / pkg) as t:
            matches = sorted(
                n
                for n in t.getnames()
                if pattern in n
                and n.endswith(".fits")
                and (field == "HD202628" or ("cont" in n and "pbcor" in n))
            )
            filename = matches[0]
            with fits.open(io.BytesIO(t.extractfile(filename).read())) as hd:
                h = hd[0].header
                img = np.squeeze(hd[0].data).astype(float) * 1e6
        w = WCS(h).celestial
        if h["BUNIT"].lower() != "jy/beam":
            raise ValueError(h["BUNIT"])
        date = Time(h["DATE-OBS"])
        beam = (h["BMAJ"] * 3600, h["BMIN"] * 3600)
        gaia = next(
            r
            for r in csv.DictReader(
                (here / "catalogs" / (field + "_GaiaDR3.csv")).open()
            )
            if r["source_id"] == GAIA[field]
        )
        host = SkyCoord(
            float(gaia["ra"]), float(gaia["dec"]), unit="deg"
        ).spherical_offsets_by(*(apparent(gaia, date.mjd) * u.arcsec))

        def measure(key, sky):
            xx, yy = w.world_to_pixel(sky)
            xx, yy = float(xx), float(yy)
            px = w.proj_plane_pixel_scales()[0].to_value(u.arcsec)
            y, x = np.indices(img.shape)
            r = np.hypot(x - xx, y - yy) * px
            local = r < np.mean(beam) / 2
            noise = img[(r >= 1.8 * np.max(beam)) & (r < 4 * np.max(beam))]
            mad = 1.4826 * np.nanmedian(abs(noise - np.nanmedian(noise)))
            flux = map_coordinates(
                img, [[yy], [xx]], order=1, mode="constant", cval=np.nan
            )[0]
            results.append(
                dict(
                    **key,
                    map_member=filename,
                    date=h["DATE-OBS"],
                    bmaj_arcsec=beam[0],
                    bmin_arcsec=beam[1],
                    ra=sky.ra.deg,
                    dec=sky.dec.deg,
                    flux_uJy_beam=flux,
                    local_MAD_uJy_beam=mad,
                    peak_half_beam_uJy=np.nanmax(img[local]) if local.any() else np.nan,
                )
            )

        for s in [r for r in cats if r["catalog"] == "GaiaDR3" and r["field"] == field]:
            ra = SkyCoord(float(s["ra"]), float(s["dec"]), unit="deg")
            peak = next(
                (
                    r
                    for r in peaks
                    if r["field"] == field
                    and r["file"] == s["file"]
                    and r["id"] == s["id"]
                ),
                None,
            )
            m = float(peak["mjd"]) if peak else 60463.09718005787
            for track, offset in [
                ("stationary", np.zeros(2)),
                ("comoving", apparent(gaia, date.mjd) - apparent(gaia, m)),
            ]:
                measure(
                    dict(field=field, file=s["file"], id=s["id"], track=track),
                    ra.spherical_offsets_by(*(offset * u.arcsec)),
                )
        if field == "HD202628":
            for r in csv.DictReader(
                (here / "catalogs" / (field + "_SIMBAD.csv")).open()
            ):
                if "FKS2019" in r["main_id"]:
                    measure(
                        dict(
                            field=field,
                            file="Faramaz2019_Table4",
                            id=r["main_id"],
                            track="stationary",
                        ),
                        SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg"),
                    )
    write_table(a.out / "ALMA_forced.csv", list(results[0]), results)


if __name__ == "__main__":
    main()
