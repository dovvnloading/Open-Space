"""Reproduce the SN 2009ib event extraction and positional controls."""
import argparse
import csv
import hashlib
import json
from io import BytesIO
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import urlencode
from urllib.request import urlopen, urlretrieve

import numpy as np
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.wcs import WCS
from scipy.spatial import cKDTree
from scipy.stats import binom

ROOT = Path(__file__).resolve().parent
EVENT_URL = ("https://cxc.cfa.harvard.edu/cdaftp/byobsid/5/16745/primary/"
             "acisf16745N002_evt2.fits.gz")
EVENT_SHA = "94791c9befbe183d47178fd032819b3f55fa65e0aa6c0301e68bfc0d95ce2ba3"
# Position from CHASE, discovery circular CBET 1902, ICRS.
RA = 15 * (4 + 17/60 + 39.92/3600)
DEC = -(62 + 46/60 + 38.7/3600)


def extract(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != EVENT_SHA:
        raise ValueError("The event file differs from the analyzed Chandra release.")
    with fits.open(path) as hdus:
        events = hdus["EVENTS"]
        hdr, data = events.header, events.data
        # WCS is encoded in the X and Y binary-table column descriptors.
        cols = [data.names.index(name) + 1 for name in ("x", "y")]
        wcs = WCS(naxis=2)
        wcs.wcs.ctype = [hdr[f"TCTYP{i}"] for i in cols]
        wcs.wcs.crval = [hdr[f"TCRVL{i}"] for i in cols]
        wcs.wcs.crpix = [hdr[f"TCRPX{i}"] for i in cols]
        wcs.wcs.cdelt = [hdr[f"TCDLT{i}"] for i in cols]
        x, y = wcs.all_world2pix(RA, DEC, 1)
        energy = data["energy"]
        band = (energy >= 500) & (energy <= 7000)
        dx = (data["x"][band].astype(float) - x) * abs(wcs.wcs.cdelt[0]) * 3600
        dy = (data["y"][band].astype(float) - y) * abs(wcs.wcs.cdelt[1]) * 3600
        return radial_counts(np.hypot(dx, dy))


EDGES = np.array([0., 1., 1.5, 2., 2.5, 3., 4., 5., 7., 8.])


def radial_counts(radius):
    return np.array([np.sum((radius >= low) & (radius < high))
                     for low, high in zip(EDGES[:-1], EDGES[1:])])


def photometry(bins):
    out = []
    for aperture, inner, outer in [(1.5, 3, 8), (1, 3, 8), (2, 3, 8),
                                    (1.5, 2.5, 5), (1.5, 4, 7)]:
        n = int(np.sum(bins[EDGES[1:] <= aperture]))
        b = int(np.sum(bins[(EDGES[:-1] >= inner) & (EDGES[1:] <= outer)]))
        ratio = aperture**2 / (outer**2 - inner**2)
        out.append(dict(radius_arcsec=aperture, annulus_arcsec=[inner, outer],
                        counts=n, background_counts=b, background_expected=b*ratio,
                        local_binomial_tail=float(binom.sf(n-1, n+b, ratio/(1+ratio)))))
    return out


def verify_catalog():
    with (ROOT / "csc-neighborhood.csv").open() as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames
        saved = {r["name"]: r for r in reader}
    q = (f"SELECT {','.join(fields)} FROM csc WHERE 1=CONTAINS("
         f"POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},0.049))")
    endpoint = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync?"
    params = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="votable", QUERY=q)
    with urlopen(endpoint + urlencode(params), timeout=90) as response:
        live = parse_single_table(BytesIO(response.read())).to_table()
    if set(live["name"]) != set(saved):
        raise ValueError("CSC neighborhood membership changed.")
    for row in live:
        original = saved[row["name"]]
        for field in fields[1:]:
            if field.endswith("_flag"):
                if str(row[field]) != original[field]:
                    raise ValueError(f"{field} changed for {row['name']}")
            else:
                old = float(original[field]) if original[field] else np.nan
                new = np.nan if np.ma.is_masked(row[field]) else float(row[field])
                np.testing.assert_allclose(new, old, atol=1e-23, rtol=1e-7)


def geometry():
    with (ROOT / "csc-neighborhood.csv").open() as fh:
        rows = list(csv.DictReader(fh))
    scale = 3600 * np.cos(np.deg2rad(DEC))
    xy = np.array([[(float(r["ra"]) - RA) * scale,
                    (float(r["dec"]) - DEC) * 3600] for r in rows])
    target = next(i for i, r in enumerate(rows) if r["name"] == "2CXO J041739.9-624638")
    offset = float(np.hypot(*xy[target]))
    # Center from Ma et al. 2023, Fig. 1. Match projected galactocentric radius.
    center = [(15*(4+17/60+35.77/3600) - RA)*scale,
              (-(62+47/60+1.2/3600) - DEC)*3600]
    rho = np.hypot(*center)
    tree = cKDTree(xy)
    controls = []
    for width in [5, 10, 20]:
        rng = np.random.default_rng(19682009)
        r = np.sqrt(rng.uniform((rho-width)**2, (rho+width)**2, 100000))
        phi = rng.uniform(0, 2*np.pi, len(r))
        points = np.column_stack((r*np.cos(phi), r*np.sin(phi))) + center
        d = tree.query(points)[0]
        controls.append(dict(radial_half_width_arcsec=width, trials=len(r),
                             within_measured_offset=int(np.sum(d <= offset)),
                             within_major_95_radius=int(np.sum(d <= float(rows[target]["error_ellipse_r0"])))))
    flux = float(rows[target]["b_flux_ap"])
    luminosity = 4*np.pi*(19.8e6*3.085677581491367e18)**2 * flux
    return dict(offset_arcsec=offset, galactocentric_radius_arcsec=float(rho),
                flux_erg_s_cm2=flux, luminosity_at_19p8_Mpc_erg_s=float(luminosity),
                controls=controls)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event", type=Path, help="Local unmodified level-2 FITS file")
    parser.add_argument("--download", action="store_true", help="Download the original event list")
    args = parser.parse_args()
    record = np.loadtxt(ROOT / "radial-counts.csv", delimiter=",", skiprows=1)
    np.testing.assert_array_equal(record[:, 0], EDGES[:-1])
    np.testing.assert_array_equal(record[:, 1], EDGES[1:])
    archived = record[:, 2].astype(int)
    if args.event or args.download:
        with TemporaryDirectory() as tmp:
            p = args.event or Path(tmp) / "events.fits.gz"
            if not args.event:
                urlretrieve(EVENT_URL, p)
            actual = extract(p)
            np.testing.assert_array_equal(actual, archived)
    if args.download:
        verify_catalog()
    print(json.dumps(dict(photometry=photometry(archived), geometry=geometry()), indent=2))


if __name__ == "__main__":
    main()
