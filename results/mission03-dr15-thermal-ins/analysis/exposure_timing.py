"""Exposure-corrected local background and four-quarter time tests for a PPS event list."""

import argparse
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.stats import binomtest, chisquare
from event_apertures import event_wcs


def aperture_ratio(path, ra, dec):
    with fits.open(path) as file:
        im = file[0].data.astype(float)
        h = file[0].header
        x, y = WCS(h, naxis=2).all_world2pix([[ra, dec]], 0)[0]
        yi, xi = np.mgrid[int(y) - 15 : int(y) + 16, int(x) - 15 : int(x) + 16]
        exposure = im[yi, xi]
        pixel_scale = abs(h["CDELT1"]) * 3600
        # Integrate the subpixel area; the 4-arcsec exposure pixels otherwise bias small apertures.
        sums = np.zeros(2)
        for dx in np.arange(-0.4375, 0.5, 0.125):
            for dy in np.arange(-0.4375, 0.5, 0.125):
                radius = np.hypot(xi + dx - x, yi + dy - y) * pixel_scale
                sums += [
                    exposure[radius < 15].sum(),
                    exposure[(radius >= 25) & (radius < 45)].sum(),
                ]
        return float(sums[0] / sums[1]) if sums[1] else None


def audit(p, ra, dec):
    e_path = p.with_name(
        p.name.replace("PIEVLI0000", "EXPMAP1000").replace("MIEVLI0000", "EXPMAP1000")
    )
    ratio = aperture_ratio(e_path, ra, dec)
    if ratio is None or ratio <= 0:
        return {"file": p.name, "covered": False}
    with fits.open(p) as f:
        hd = f["EVENTS"]
        e = hd.data
        h = hd.header
        x, y = event_wcs(h).all_world2pix([[ra, dec]], 1)[0]
        r = np.hypot(
            (e["X"] - x) * h["REFXCDLT"] * 3600, (e["Y"] - y) * h["REFYCDLT"] * 3600
        )
        pn = h["INSTRUME"] == "EPN"
        clean = (
            (e["FLAG"] == 0)
            & (e["PATTERN"] <= (4 if pn else 12))
            & (e["PI"] >= 200)
            & (e["PI"] < 1000)
        )
        core = clean & (r < 15)
        if not core.any():
            return {"file": p.name, "covered": True, "count": 0}
        # Identify the CCD containing the source by the mode among aperture events.
        ccd = np.bincount(e["CCDNR"][core]).argmax()
        source_gtis = f[f"STDGTI{ccd:02}"].data
        t = e["TIME"]
        live = np.zeros(len(t))
        in_gti = np.zeros(len(t), dtype=bool)
        base = 0.0
        for g in source_gtis:
            m = (t >= g["START"]) & (t <= g["STOP"])
            live[m] = base + t[m] - g["START"]
            in_gti |= m
            base += g["STOP"] - g["START"]
        if not np.all(e["CCDNR"][core] == ccd):
            raise ValueError("Time audit assumes the source occupies a single CCD")
        n = np.sum(core & in_gti)
        nback = np.sum(clean & (r >= 25) & (r < 45) & in_gti)
        time_counts = np.histogram(live[core & in_gti], bins=np.linspace(0, base, 5))[0]
        return {
            "file": p.name,
            "covered": True,
            "energy_keV": [0.2, 1.0],
            "source_ccd": int(ccd),
            "ccd_gti_seconds": base,
            "aperture_counts": int(n),
            "annulus_counts": int(nback),
            "exposure_area_ratio": ratio,
            "estimated_background": ratio * nback,
            "conditional_p_equal_surface_rate": float(
                binomtest(
                    int(n), int(n + nback), ratio / (1 + ratio), alternative="greater"
                ).pvalue
            ),
            "four_equal_gti_quarter_counts": time_counts.tolist(),
            "constant_total_rate_chi2_p": float(chisquare(time_counts).pvalue),
        }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("ra", type=float)
    p.add_argument("dec", type=float)
    p.add_argument("directory", type=Path)
    a = p.parse_args()
    print(
        json.dumps(
            [
                audit(f, a.ra, a.dec)
                for f in sorted(a.directory.glob("*[MP]IEVLI*.FTZ"))
            ],
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
