"""Replay the published image-center audit with the Python standard library."""

import csv
import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = math.pi / 180


def angsep(lon, lat, lon2, lat2):
    a = math.sin((lat2 - lat) * D / 2) ** 2
    a += math.cos(lat * D) * math.cos(lat2 * D) * math.sin((lon2 - lon) * D / 2) ** 2
    return 2 * math.asin(math.sqrt(min(1, a))) / D


def native_wcs(x, y, f):
    # FITS SIN inverse. The release uses diagonal axes, zero rotation/slant.
    lng0 = f("CRVAL1") * D
    lat0 = f("CRVAL2") * D
    u = (x - f("CRPIX1")) * f("CDELT1") * D
    v = (y - f("CRPIX2")) * f("CDELT2") * D
    z = math.sqrt(1 - u * u - v * v)
    lng = lng0 + math.atan2(u, z * math.cos(lat0) - v * math.sin(lat0))
    lat = math.asin(z * math.sin(lat0) + v * math.cos(lat0))
    return lng / D % 360, lat / D


def pixel(lon, lat, f):
    # Direct SIN projection. FITS pixel centers are one-based.
    l = (lon - f("CRVAL1")) * D
    b = lat * D
    b0 = f("CRVAL2") * D
    u = math.cos(b) * math.sin(l)
    v = math.sin(b) * math.cos(b0) - math.cos(b) * math.sin(b0) * math.cos(l)
    return f("CRPIX1") + u / D / f("CDELT1"), f("CRPIX2") + v / D / f("CDELT2")


def main():
    rows = []
    for row in csv.DictReader((HERE / "headers.csv").open()):

        def f(key, row=row):
            return float(row[key])

        if (
            (row["CTYPE1"], row["CTYPE2"]) != ("GLON-SIN", "GLAT-SIN")
            or f("CROTA1")
            or f("CROTA2")
        ):
            raise ValueError("Unexpected coordinate system for " + row["filename"])
        nx, ny = f("nx"), f("ny")
        if nx != f("NAXIS1") or ny != f("NAXIS2"):
            raise ValueError("Differing image sizes")
        x, y = (nx + 1) / 2, (ny + 1) / 2
        lon, lat = native_wcs(x, y, f)
        lx, ly = pixel(f("listed_l"), f("listed_b"), f)
        sep = angsep(lon, lat, f("listed_l"), f("listed_b")) * 60
        # Empirical model evaluated separately from the exact native WCS.
        ml = (f("CRVAL1") - (x - f("CRPIX1")) * f("CDELT1")) % 360
        mb = f("CRVAL2") + (y - f("CRPIX2")) * f("CDELT2")
        residual = 3600 * angsep(ml, mb, f("listed_l"), f("listed_b"))
        inside = 0.5 <= lx <= nx + 0.5 and 0.5 <= ly <= ny + 0.5
        rows.append(
            {
                "filename": row["filename"],
                "listed_l": row["listed_l"],
                "listed_b": row["listed_b"],
                "native_l": round(lon, 8),
                "native_b": round(lat, 8),
                "offset_arcmin": round(sep, 5),
                "listed_x_fits": round(lx, 3),
                "listed_y_fits": round(ly, 3),
                "listed_inside": int(inside),
                "model_residual_arcsec": round(residual, 5),
            }
        )
    with (HERE / "coordinates.csv").open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    seps = [r["offset_arcmin"] for r in rows]
    residuals = [r["model_residual_arcsec"] for r in rows]
    print(
        json.dumps(
            {
                "images": len(rows),
                "outside": sum(1 - r["listed_inside"] for r in rows),
                "within_5_arcmin": sum(s <= 5 for s in seps),
                "median_arcmin": statistics.median(seps),
                "max_arcmin": max(seps),
                "min_arcmin": min(seps),
                "median_model_residual_arcsec": statistics.median(residuals),
                "max_model_residual_arcsec": max(residuals),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
