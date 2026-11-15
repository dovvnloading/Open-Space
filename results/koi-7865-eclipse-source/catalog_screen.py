"""Reproduce the KOI to Gaia EB screen. Run from this directory with Python 3.10+."""
import hashlib
import io
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
CACHE = HERE / ".cache"
KOI = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urlencode({
    "query": "select kepid,kepoi_name,kepler_name,koi_period,koi_period_err1,koi_time0bk,koi_depth,koi_duration,koi_disposition,koi_pdisposition,koi_fpflag_co,koi_fpflag_ec,koi_fpflag_ss,koi_fpflag_nt,koi_kepmag,ra,dec from cumulative",
    "format": "csv",
})
VIZ = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
EB = VIZ + urlencode({
    "-source": "I/358/veb", "-c": "292.25 44.5", "-c.rd": "12",
    "-out.max": "unlimited", "-out": "Source,RA_ICRS,DE_ICRS,Freq",
})
NEARBY = VIZ + urlencode({
    "-source": "I/355/gaiadr3", "-c": "296.18478 43.980751", "-c.rs": "20",
    "-out": "Source,RA_ICRS,DE_ICRS,Gmag,pmRA,pmDE", "-out.max": "1000",
})
HARMONICS = np.array([0.5, 1., 2., 3., 4., 1./3, 0.25])


def fetch(name, url, sha256=None):
    CACHE.mkdir(exist_ok=True)
    path = CACHE / name
    if not path.exists():
        data = urlopen(url, timeout=120).read()
        if sha256 and hashlib.sha256(data).hexdigest() != sha256:
            raise ValueError(f"Unexpected download checksum for {name}")
        if len(data) < 100 or data.lstrip().lower().startswith(b'<!doctype'):
            raise ValueError(f"Unexpected download body for {name}")
        path.write_bytes(data)
    if sha256 and hashlib.sha256(path.read_bytes()).hexdigest() != sha256:
        raise ValueError(f"Unexpected checksum for {path}")
    return path


def read_vizier(path):
    rows = pd.read_csv(path, sep="\t", comment="#", dtype=str).iloc[2:].reset_index(drop=True)
    if len(rows) == 0:
        raise ValueError(f"No catalog rows in {path}")
    return rows


def unit_vectors(ra, dec):
    lon, lat = np.deg2rad(np.asarray(ra, float)), np.deg2rad(np.asarray(dec, float))
    return np.column_stack((np.cos(lon)*np.cos(lat), np.sin(lon)*np.cos(lat), np.sin(lat)))


def search(koi, eb, shift_arcsec=0):
    # Offsets shift declination, preserving periods, angular clustering and magnitude.
    stars = unit_vectors(koi.ra, koi.dec + shift_arcsec/3600)
    eb_xyz = unit_vectors(eb.RA_ICRS, eb.DE_ICRS)
    tree = cKDTree(eb_xyz)
    frequencies = eb.Freq.astype(float).to_numpy()
    matches = []
    for i, neighbors in enumerate(tree.query_ball_point(stars, 2*np.sin(np.deg2rad(90/3600)/2))):
        candidate = koi.iloc[i]
        for j in neighbors:
            ratio = candidate.koi_period*frequencies[j]
            harmonic = HARMONICS[np.argmin(abs(ratio/HARMONICS-1))]
            mismatch = ratio/harmonic - 1
            if abs(mismatch) > .002:
                continue
            sep = np.rad2deg(2*np.arcsin(np.linalg.norm(stars[i]-eb_xyz[j])/2))*3600
            matches.append(dict(koi=candidate.kepoi_name, disposition=candidate.koi_disposition,
                                gaia_id=eb.Source.iloc[j], separation_arcsec=sep,
                                period_ratio_harmonic=harmonic, fractional_mismatch=mismatch,
                                koi_period_d=candidate.koi_period, gaia_period_d=1/frequencies[j],
                                koi_depth_ppm=candidate.koi_depth, stellar_eclipse_flag=candidate.koi_fpflag_ss,
                                centroid_flag=candidate.koi_fpflag_co, ephemeris_flag=candidate.koi_fpflag_ec))
    return pd.DataFrame(matches)


def main():
    koi = pd.read_csv(fetch("koi.csv", KOI))
    eb = read_vizier(fetch("eb.tsv", EB))
    matches = search(koi, eb)
    active = matches[matches.disposition != "FALSE POSITIVE"].sort_values("separation_arcsec")
    active.to_csv(HERE/"active_matches.csv", index=False)
    radius = 6.
    tolerance = 1e-5
    keep = (active.separation_arcsec <= radius) & (abs(active.fractional_mismatch) < tolerance)
    print(f"{len(koi)} KOIs and {len(eb)} Gaia EBs; {len(matches)} broad matches, {len(active)} active; {keep.sum()} within {radius} arcsec and {tolerance} fractional period mismatch.")
    # Descriptive controls. The tight cut summarizes the identified lead.
    counts = []
    for shift in [-3600, -3000, -2400, -1800, -1200, -600, 600, 1200, 1800, 2400, 3000, 3600]:
        control = search(koi[koi.koi_disposition != "FALSE POSITIVE"], eb, shift)
        tight = ((control.separation_arcsec <= radius) & (abs(control.fractional_mismatch) < tolerance)).sum() if len(control) else 0
        counts.append({"declination_shift_arcsec": shift, "broad_matches": len(control), "tight_matches": tight})
    pd.DataFrame(counts).to_csv(HERE/"offset_controls.csv", index=False)
    print("Shifted wide/tight counts:", sum(i["broad_matches"] for i in counts), sum(i["tight_matches"] for i in counts))


if __name__ == "__main__":
    main()
