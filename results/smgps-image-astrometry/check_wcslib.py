"""Independently verify the standard-library projection using WCSLIB via Astropy."""

import csv

from astropy.io.fits import Header
from astropy.wcs import WCS
from audit import HERE, angsep


def main():
    with (HERE / "coordinates.csv").open() as stream:
        results = {r["filename"]: r for r in csv.DictReader(stream)}
    maxerr = 0
    with (HERE / "headers.csv").open() as stream:
        for r in csv.DictReader(stream):
            h = Header()
            h["NAXIS"] = 2
            for k in [
                "NAXIS1",
                "NAXIS2",
                "CRPIX1",
                "CRPIX2",
                "CRVAL1",
                "CRVAL2",
                "CDELT1",
                "CDELT2",
                "CTYPE1",
                "CTYPE2",
                "CROTA1",
                "CROTA2",
            ]:
                if k.startswith("NAXIS"):
                    h[k] = int(r[k])
                else:
                    h[k] = r[k] if k.startswith("CTYPE") else float(r[k])
            w = WCS(h, naxis=2)
            lon, lat = w.pixel_to_world_values(
                (h["NAXIS1"] - 1) / 2, (h["NAXIS2"] - 1) / 2
            )
            out = results[r["filename"]]
            maxerr = max(
                maxerr,
                angsep(lon, lat, float(out["native_l"]), float(out["native_b"])) * 3600,
            )
    if maxerr > 0.001:
        raise ValueError(f"Centers differ by {maxerr} arcsec")
    print(
        f"{len(results)} centers verified via WCSLIB. Maximum separation {maxerr:.6f} arcsec."
    )


if __name__ == "__main__":
    main()
