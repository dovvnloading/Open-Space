"""Check local translational offsets of the native DR15 coordinates using Gaia stars."""

import argparse
import csv
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from scipy.spatial import cKDTree


def measure(path, target, cache):
    with fits.open(path, memmap=True) as h:
        data = h[1].data
        row = np.where(data["SRCID"] == target)[0].item()
        center = SkyCoord(data["RA"][row], data["DEC"][row], unit="deg")
        epoch = Time(float(data["MJD_FIRST"][row]), format="mjd").jyear
        radius = 900.0
        rr, dd = center.ra.deg, center.dec.deg
        box = (
            abs((data["RA"] - rr + 180) % 360 - 180) * np.cos(np.deg2rad(dd)) < 0.26
        ) & (abs(data["DEC"] - dd) < 0.26)
        # Only use objects whose stack was assigned the target's observation prefix.
        idx = np.where(
            box
            & (data["SRCID"] // 100000 == target // 100000)
            & (data["SRCID"] != target)
            & (data["SUM_FLAG"] <= 1)
            & (data["EXTENT"] == 0)
            & (data["DIST_NN"] >= 12)
            & (data["RADEC_ERR"] < 2.0)
            & (data["STACK_DET_ML"] >= 20)
        )[0]
        x = SkyCoord(data["RA"][idx], data["DEC"][idx], unit="deg")
        keep = center.separation(x).arcsec < radius
        x = x[keep]
        ids = data["SRCID"][idx[keep]]
    fields = "source_id,ra,dec,pmra,pmdec,phot_g_mean_mag,parallax,parallax_over_error"
    query = f"SELECT {fields} FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{rr},{dd},0.255)) AND phot_g_mean_mag < 21"
    url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urlencode(
        dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query)
    )
    raw = cache / f"{target}-gaia-field.csv"
    if not raw.exists():
        with urlopen(url, timeout=90) as r:
            raw.write_bytes(r.read())
    cat = list(csv.DictReader(io.StringIO(raw.read_text())))
    if not cat:
        raise ValueError("Missing Gaia field")

    def values(name, default=np.nan):
        return np.array([float(z.get(name) or default) for z in cat])

    dec = values("dec") + values("pmdec", 0) * (epoch - 2016) / 3.6e6
    ra = values("ra") + values("pmra", 0) * (epoch - 2016) / 3.6e6 / np.cos(
        np.deg2rad(dec)
    )
    g = SkyCoord(ra, dec, unit="deg")
    frame = center.skyoffset_frame()

    def coords(c):
        d = c.transform_to(frame)
        return np.column_stack((d.lon.arcsec, d.lat.arcsec))

    gp, xp = coords(g), coords(x)
    tree = cKDTree(gp)
    search, neigh = tree.query(xp, k=2)
    pairs = []
    for i, (dist, ii) in enumerate(zip(search, neigh)):
        j = ii[0]
        if dist[0] > 3.5 or dist[1] < 7 or cKDTree(xp).query(gp[j])[1] != i:
            continue
        pairs.append(
            dict(
                srcid=str(ids[i]),
                gaia_id=cat[j]["source_id"],
                g=float(cat[j]["phot_g_mean_mag"]),
                sep=round(dist[0], 3),
                dx=float(gp[j, 0] - xp[i, 0]),
                dy=float(gp[j, 1] - xp[i, 1]),
            )
        )
    # Empirical control for false matches. Displace every X-ray position by 30 arcsec.
    ncontrol = []
    for angle in np.arange(8) * np.pi / 4:
        z, _ = tree.query(xp + 30 * np.array([np.cos(angle), np.sin(angle)]), k=2)
        ncontrol.append(int(np.sum((z[:, 0] <= 3.5) & (z[:, 1] >= 7))))
    res = dict(
        target=str(target),
        epoch=epoch,
        gaia_query=url,
        gaia_in_field=len(cat),
        dr15_references_tested=len(x),
        pairs=pairs,
        shifted_control_counts=ncontrol,
    )
    if pairs:
        res.update(
            median_east_arcsec=float(np.median([t["dx"] for t in pairs])),
            median_north_arcsec=float(np.median([t["dy"] for t in pairs])),
        )
    return res


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fits")
    p.add_argument("srcid", type=int)
    p.add_argument("cache", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    a.cache.mkdir(exist_ok=True, parents=True)
    z = measure(a.fits, a.srcid, a.cache)
    a.output.write_text(json.dumps(z, indent=2) + "\n")
    print(z)
