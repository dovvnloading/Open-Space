"""Stack per-eclipse difference images and locate each clock relative to Gaia."""
from pathlib import Path
import json
import sys
import io
from zipfile import ZipFile
from urllib.request import urlopen

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import pandas as pd

from extract import digest
from fit_time import CLOCKS, HERE

HOST = (126.771492341, -50.765804482)
NEIGHBOR = (126.765395868, -50.767146484)


def centroid(image, wcs, hw):
    x, y = wcs.world_to_pixel_values(*HOST)
    yy, xx = np.indices(image.shape)
    rim = (abs(xx-x) > 2.1) | (abs(yy-y) > 2.1)
    weights = np.maximum(image - np.median(image[rim]), 0)
    weights *= (abs(xx-x) <= hw) & (abs(yy-y) <= hw)
    if not weights.sum():
        return np.nan, np.nan
    ra, dec = wcs.pixel_to_world_values(np.sum(xx*weights)/weights.sum(),
                                         np.sum(yy*weights)/weights.sum())
    return float((ra-HOST[0])*np.cos(np.deg2rad(HOST[1]))*3600), float((dec-HOST[1])*3600)


def measure(hdus, sector, visits):
    d = hdus[1].data
    keep = (d["QUALITY"] == 0) & np.isfinite(d["TIME"]) & np.isfinite(d["FLUX"]).all(axis=(1, 2))
    cube = d["FLUX"][keep].astype(float)
    t = pd.read_csv(HERE / "data" / f"s{sector}.csv.gz").time.to_numpy()
    if len(t) != len(cube):
        raise ValueError("Cadences in cube and barycentric time series differ")
    wcs = WCS(hdus[2].header)
    sector_visits = visits[visits.sector == sector]
    ephemerides = {r.clock: (CLOCKS[r.clock][0], CLOCKS[r.clock][1] + r.dt/1440)
                   for r in sector_visits.itertuples()}
    for name, (period, epoch) in ephemerides.items():
        keep = np.ones(len(t), dtype=bool)
        for label, (other_p, other_e) in ephemerides.items():
            if label != name:
                other_phase = ((t-other_e)/other_p+.5) % 1-.5
                keep &= abs(other_phase) > 4 * CLOCKS[label][2]
        images = []
        cycles = np.arange(np.floor((t.min()-epoch)/period), np.ceil((t.max()-epoch)/period)+1)
        for cycle in cycles:
            delta = t - (epoch + period*cycle)
            inside = keep & (abs(delta) < .035)
            before = keep & (delta > -.24) & (delta < -.07)
            after = keep & (delta < .24) & (delta > .07)
            if min(inside.sum(), before.sum(), after.sum()) >= 2:
                baseline = (np.median(cube[before],axis=0) + np.median(cube[after],axis=0))/2
                images.append(baseline - np.median(cube[inside],axis=0))
        if not images:
            continue
        image = np.median(images,axis=0)
        for size in [1.8, 2.3]:
            east, north = centroid(image, wcs, size)
            yield dict(clock=name, sector=sector, eclipses=len(images), radius_px=size,
                       east_arcsec=east, north_arcsec=north, brightest_pixel_eps=image.max())


if __name__ == "__main__":
    visits = pd.read_csv(HERE / "visits.csv")
    visits = visits[(visits.aperture == "box") & (visits.window == .17)]
    manifest = json.loads((HERE / "sources.json").read_text())
    rows = []
    for row in manifest["tess"]:
        s = row["sector"]
        if s not in set(visits.sector):
            continue
        if len(sys.argv) > 1:
            path = Path(sys.argv[1]) / f"s{s}.fits"
        else:
            raw = urlopen(row["cutout"], timeout=240).read()
            if raw.startswith(b"PK"):
                with ZipFile(io.BytesIO(raw)) as archive:
                    member, = (n for n in archive.namelist() if n.endswith(".fits"))
                    raw = archive.read(member)
            path = io.BytesIO(raw)
        with fits.open(path) as hdus:
            if digest(hdus) != row["arrays_sha256"]:
                raise ValueError(f"Cutout {s} changed")
            rows.extend(measure(hdus,s,visits))
    pd.DataFrame(rows).to_csv(HERE / "positions.csv", index=False, float_format="%.6f")
