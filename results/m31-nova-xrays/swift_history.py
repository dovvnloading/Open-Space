"""Measure the Swift PC energy distribution toward M31N 2021-06b over 58 visits."""
import argparse
import csv
import sys

import numpy as np
from scipy.stats import binom

from xray_bands import product, read
from pathlib import Path

# Swift PC calibrated PI bins (10 eV/channel), with inclusive 10 keV endpoint.
BANDS = ((30, 40), (40, 50), (50, 60), (60, 70), (70, 85), (85, 100),
         (30, 100), (100, 200), (200, 1001))
FIELDS = ("obsid", "date_utc", "period", "center", "radius_pix", "pi_min", "pi_max_inclusive",
          "aperture_exposure_s", "on", "off", "alpha", "east_sum_arcsec", "north_sum_arcsec")


def measure(cache):
    from astropy.io import fits
    from astropy.wcs import WCS

    cache.mkdir(parents=True, exist_ok=True)
    nova = next(t for t in read("targets.csv") if t["name"] == "M31N 2021-06b")
    neighbor = next(r for r in read("csc-neighbors.csv") if r["name"] == "2CXO J004258.0+411809")
    coords = [[float(r[a]), float(r[b])] for r, a, b in
              ((nova, "ra_deg", "dec_deg"), (neighbor, "ra", "dec"))]
    for obs in read("swift-timeline-obs.csv"):
        evt = product(cache, obs["event_url"], obs["event_sha256"])
        exp = product(cache, obs["map_url"], obs["map_sha256"])
        with fits.open(evt) as hd, fits.open(exp) as maps:
            events = hd["EVENTS"]
            image = maps[0].data
            header = maps[0].header
            for axis, colname in enumerate(("X", "Y"), 1):
                col = events.columns.names.index(colname) + 1
                for ek, mk in (("TCRPX", "CRPIX"), ("TCRVL", "CRVAL"), ("TCDLT", "CDELT")):
                    if not np.isclose(events.header[f"{ek}{col}"], header[f"{mk}{axis}"]):
                        raise ValueError("Event/map WCS mismatch")
            (x, y), (nx, ny) = WCS(header).all_world2pix(coords, 1)
            date = events.header["DATE-OBS"][:10]
            period = "pre" if date < "2021-06-12" else "summer" if date < "2021-10-01" else "late"
            # Cover both targets and the ring of control apertures with a margin.
            x0, y0 = max(1, int(x)-43), max(1, int(y)-43)
            im = image[y0-1:int(y)+44, x0-1:int(x)+44]
            gy, gx = np.ogrid[y0:y0+im.shape[0], x0:x0+im.shape[1]]
            d = events.data
            good = (d["GRADE"] >= 0) & (d["GRADE"] <= 12)
            centers = [("nova", x, y), ("neighbor", nx, ny)]
            for i in range(8):
                angle = i * np.pi / 4
                centers.append((f"shift_{i}", x + 12*np.cos(angle), y + 12*np.sin(angle)))
            for name, cx, cy in centers:
                dr = np.hypot(d["X"]-cx, d["Y"]-cy)
                pixels = np.hypot(gx-cx, gy-cy)
                sky = (pixels > 12) & (pixels <= 24)
                for radius in ((4, 6, 8) if name == "nova" else (6,)):
                    source = pixels <= radius
                    time = im[source].mean(dtype=float)
                    if time < 30 or im[sky].mean(dtype=float) < 30:
                        continue
                    alpha = im[source].sum(dtype=float) / im[sky].sum(dtype=float)
                    bands = BANDS if name == "nova" and radius == 6 else ((30, 100),)
                    for low, high in bands:
                        selected = good & (d["PI"] >= low) & (d["PI"] < high)
                        on = selected & (dr <= radius)
                        off = selected & (dr > 12) & (dr <= 24)
                        east = (d["X"][on]-x).sum() * header["CDELT1"] * 3600
                        north = (d["Y"][on]-y).sum() * header["CDELT2"] * 3600
                        yield (obs["obsid"], date, period, name, radius, low, high-1,
                               f"{time:.6f}", on.sum(), off.sum(), f"{alpha:.9f}",
                               f"{east:.5f}", f"{north:.5f}")


def summarize():
    rows = read("swift-timeline-counts.csv")
    print("period center radius PI Nobs exposure_ks on scaled_background east_mean north_mean")
    for period in ("pre", "summer", "late"):
        groups = {}
        for r in rows:
            if r["period"] != period:
                continue
            key = (r["center"], int(r["radius_pix"]), int(r["pi_min"]), int(r["pi_max_inclusive"]))
            groups.setdefault(key, []).append(r)
        for (cen, rad, lo, hi), group in sorted(groups.items()):
            on = sum(int(r["on"]) for r in group)
            bg = sum(int(r["off"]) * float(r["alpha"]) for r in group)
            time = sum(float(r["aperture_exposure_s"]) for r in group)
            xy = [sum(float(r[k]) for r in group)/on if on else np.nan
                  for k in ("east_sum_arcsec", "north_sum_arcsec")]
            print(f"{period} {cen} {rad} {lo}-{hi} {len(group)} {time/1000:.3f} "
                  f"{on} {bg:.2f} {xy[0]:.2f} {xy[1]:.2f}")
    terms = [r for r in rows if r["center"] == "nova" and int(r["radius_pix"]) == 6
             and int(r["pi_min"]) == 30 and int(r["pi_max_inclusive"]) == 99]
    for r in terms:
        on, off, alpha = int(r["on"]), int(r["off"]), float(r["alpha"])
        p = binom.sf(on-1, on+off, alpha/(1+alpha))
        print(f"{r['date_utc']} {r['obsid']} soft on={on} off={off} expected={alpha*off:.2f} local_p={p:.3g}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--extract", type=Path, metavar="CACHE", help="Download and count archived FITS")
    args = p.parse_args()
    if args.extract:
        writer = csv.writer(sys.stdout, lineterminator="\n")
        writer.writerow(FIELDS)
        writer.writerows(measure(args.extract))
    else:
        summarize()
