"""Rebuild the WISE epochs and local frame from the public archive products."""
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import argparse
import hashlib

import numpy as np
import pandas as pd
from astropy.coordinates import get_body_barycentric
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
RA, DEC = 300.209094, 36.497796
EAST_SCALE = np.cos(np.deg2rad(DEC)) * 3.6e6
ADQL = ("SELECT source_id,ra,dec,pmra,pmdec,parallax,phot_g_mean_mag,bp_rp,ruwe FROM gaiadr3.gaia_source "
        "WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',300.209094,36.49805,0.04)) "
        "AND phot_g_mean_mag<18.5")


def cached(cache, name, url):
    path = cache / name
    if not path.exists():
        with urlopen(url, timeout=180) as response:
            buffer = response.read()
        temporary = path.with_suffix(".part")
        temporary.write_bytes(buffer)
        temporary.replace(path)
    return path


def measure_product(product, data, refs, output):
    a = data.sort_values("mjd").reset_index(drop=True).copy()
    if product == "L1b":
        a["epoch"] = (a.mjd.diff() > 30).cumsum()
    a["x"] = (a.ra-RA)*EAST_SCALE
    a["y"] = (a.dec-DEC)*3.6e6
    t = Time(a.mjd, format="mjd")
    a["year"] = t.jyear
    ra, dec = np.deg2rad([RA, DEC])
    X, Y, Z = get_body_barycentric("earth", t).xyz.value
    a["fx"] = np.sin(ra)*X - np.cos(ra)*Y
    a["fy"] = np.cos(ra)*np.sin(dec)*X + np.sin(ra)*np.sin(dec)*Y - np.cos(dec)*Z
    # This aperture (1.5 or 1.8 arcsec) identifies the moving W2 source.
    # Centroid residuals inside it remain free in both dimensions.
    residual = np.hypot(a.x-133.4*a.fx-6.1*(a.year-2015.64), a.y-133.4*a.fy-372.8*(a.year-2015.64))
    if product == "L1b":
        lens = a[(residual < 1500) & a.w2mpro.between(12.15, 13.25)]
        radius, min_visits = 1500, 18
    else:
        lens = a[residual < 1800]
        radius, min_visits = 1200, 14
    lens_group = lens.groupby("epoch")
    visits = lens_group[["year","x","y","fx","fy"]].median()
    visits["n"] = lens_group.size()
    if product == "W2" and not (visits.n == 1).all():
        raise ValueError("Ambiguous W2 source selection")
    visits.insert(0,"product",product)
    visits.reset_index()[["product","epoch","year","x","y","fx","fy","n"]].to_csv(
        output / f"{product}-epochs.csv", index=False, float_format="%.10g")

    sep, ids = cKDTree(refs[["x", "y"]]).query(a[["x", "y"]])
    keep = sep < radius
    if product == "W2":
        keep &= a.fracflux > .5
    detections = a.loc[keep].copy()
    catalog = refs.iloc[ids[keep]].set_axis(detections.index)
    detections["source_id"] = catalog.source_id
    detections["dx"] = detections.x-catalog.x-(detections.year-2016)*catalog.pmra-detections.fx*catalog.parallax
    detections["dy"] = detections.y-catalog.y-(detections.year-2016)*catalog.pmdec-detections.fy*catalog.parallax
    offsets = detections.groupby(["source_id","epoch"])[["dx","dy"]].median().reset_index()
    offsets = offsets[offsets.groupby("source_id").epoch.transform("size") >= min_visits].copy()
    offsets.insert(0,"product",product)
    return offsets


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    params = urlencode(dict(QUERY=ADQL, FORMAT="csv", LANG="ADQL", REQUEST="doQuery"))
    gaia = cached(args.cache, "gaiafield.csv", "https://gea.esac.esa.int/tap-server/tap/sync?"+params)
    c = pd.read_csv(gaia, dtype={"source_id":str})
    c["x"] = (c.ra-RA)*EAST_SCALE
    c["y"] = (c.dec-DEC)*3.6e6
    neighbor, _ = cKDTree(c[["x","y"]]).query(c[["x","y"]], k=2)
    refs = c[(neighbor[:,1]>6500) & (c.ruwe<1.4) & (c.phot_g_mean_mag<16.5) & c.parallax.notna()]
    fields = "ra,dec,mjd,w2mpro,w2snr,w2rchi2,qual_frame,cc_flags,scan_id,frame_num,nb"
    params = urlencode(dict(catalog="neowiser_p1bs_psd", spatial="Cone", objstr="300.209094 36.49805",
                            radius=120, radunits="arcsec", outfmt=1, selcols=fields, constraints="w2snr>5"))
    field = cached(args.cache, "field.tbl", "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?"+params)
    d = Table.read(field, format="ascii.ipac").to_pandas()
    selected = d[(d.qual_frame>0)&(d.w2snr>=8)&(d.w2rchi2<2)&(d.nb==1)&(d.cc_flags.str[1]=="0")]
    offsets = [measure_product("L1b", selected, refs, args.output)]
    epochs = []
    for row in pd.read_csv(HERE/"unTimely-inputs.csv").itertuples(index=False):
        filename = cached(args.cache, f"w2-{row.epoch:02d}.fits.gz", row.url)
        if hashlib.sha256(filename.read_bytes()).hexdigest() != row.sha256:
            raise ValueError(f"unTimely checksum differs for epoch {row.epoch}")
        t = fits.getdata(filename)
        cone = np.hypot((t["ra"]-RA)*.804,t["dec"]-DEC)<.04
        good = cone & (t["qf"]>.8) & (t["flux"]>10*t["dflux"]) & (t["flags_unwise"]==0)
        sub = t[good]
        epochs.append(pd.DataFrame(dict(ra=sub["ra"].astype(float), dec=sub["dec"].astype(float),
                                        fracflux=sub["fracflux"].astype(float), mjd=row.mjd, epoch=row.epoch)))
    offsets.append(measure_product("W2", pd.concat(epochs), refs, args.output))
    pd.concat(offsets).sort_values(["product","source_id","epoch"]).to_csv(
        args.output/"reference-offsets.csv", index=False, float_format="%.7f")


if __name__ == "__main__":
    main()
