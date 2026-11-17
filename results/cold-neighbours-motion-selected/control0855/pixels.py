"""Fetch sampled W2 L1b cutouts and replay five independent fixed-sky image stacks."""
from concurrent.futures import ThreadPoolExecutor
import argparse
import csv
import hashlib
from pathlib import Path
from urllib.request import urlopen

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from native import cutout

BASE = Path(__file__).resolve().parent
FRAME_URL = "https://irsa.ipac.caltech.edu/ibe/search/wise/neowiser/p1bm_frm?POS=133.78,-7.24&SIZE=0.002"
RA, DEC = 133.773, -7.243
COSDEC = np.cos(np.deg2rad(DEC))
CHOICES = [0, 5, 10, 15, 20]
SCALE = 1.375


def make_manifest():
    table = urlopen(FRAME_URL, timeout=90).read()
    data = Table.read(table.decode().splitlines(), format="ascii.ipac").to_pandas()
    data = data[(data.band == 2) & (data.qual_frame == 10)].sort_values("mjd_obs").copy()
    data["visit"] = (data.mjd_obs.diff() > 20).cumsum()
    data["offset"] = np.hypot((data.crval1 - RA)*COSDEC, data.crval2 - DEC)
    rows = []
    for v in CHOICES:
        for f in data[data.visit == v].nsmallest(6, "offset").sort_values("mjd_obs").itertuples():
            scan, frame = str(f.scan_id), f"{f.frame_num:03d}"
            name = f"{scan}{frame}-w2-int-1b.fits"
            url = ("https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/neowiser/images/"
                   f"p1bm_frm/{scan[-2:]}/{scan}/{frame}/{name}")
            rows.append(dict(visit=v, mjd=f.mjd_obs, file=name, url=url))
    return rows


def load_image(row, cache):
    path = cache / (row["file"] + ".gz")
    if not path.exists():
        cutout(row["url"], RA, DEC).writeto(path)
    with fits.open(path) as hdus:
        header = hdus[0].header.copy()
        image = np.asarray(hdus[0].data, float)
    if abs(header["MJD_OBS"] - float(row["mjd"])) > 1.e-5:
        raise ValueError(f"Wrong timestamp {path}")
    # The two source byte-range digests retain the native FITS provenance.
    digest = hashlib.sha256(np.asarray(image, dtype=">f4").tobytes()).hexdigest()
    for field, key in [("input_fits_header_sha256","INHDRSHA"), ("input_fits_row_strip_sha256","INRAWSHA")]:
        if row.get(field) and row[field] != header[key]:
            raise ValueError(f"Changed archive input: {path}")
    if row.get("science_pixels_sha256") and row["science_pixels_sha256"] != digest:
        raise ValueError(f"Changed cutout science array: {path}")
    row["science_pixels_sha256"] = digest
    row["input_fits_header_sha256"] = header["INHDRSHA"]
    row["input_fits_row_strip_sha256"] = header["INRAWSHA"]
    return row, header, image


def mad(x):
    return np.nanmedian(np.abs(x - np.nanmedian(x)))*1.4826


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    args.cache.mkdir(exist_ok=True, parents=True)
    manifest = BASE/"image-manifest.csv"
    rows = list(csv.DictReader(manifest.open())) if manifest.exists() else make_manifest()
    with ThreadPoolExecutor(6) as pool:
        loaded = list(pool.map(lambda r: load_image(r, args.cache), rows))
    if not manifest.exists():
        pd.DataFrame(rows).to_csv(manifest, index=False)
    ep = pd.read_csv(BASE/"epochs.csv").set_index("visit")
    # East is left, north up. The same tangent-plane approximation is used by the astrometric fit.
    x = np.arange(42,-37,-1)*SCALE
    y = np.arange(-18,19)*SCALE
    xx, yy = np.meshgrid(x, y)
    ra, dec = RA+xx/(3600*COSDEC), DEC+yy/3600
    stacks, splits = {}, {}
    for v in CHOICES:
        projected = []
        for row, header, array in loaded:
            if int(row["visit"]) != v:
                continue
            ux, uy = WCS(header).celestial.world_to_pixel_values(ra, dec)
            cut = map_coordinates(array, [uy, ux], order=1, cval=np.nan)
            # Common approximate flux scale of nominal W2 zero point 20.
            projected.append((cut-np.nanmedian(cut))*10**(0.4*(20.-header["MAGZP"])))
        splits[v] = [np.nanmedian(projected[i::2], axis=0) for i in (0,1)]
        stacks[v] = np.nanmedian(projected, axis=0)
    results = []
    panels = []
    for v in CHOICES:
        background = np.nanmedian([stacks[w] for w in CHOICES if v != w], axis=0)
        dist = np.hypot(xx-ep.loc[v,"x"], yy-ep.loc[v,"y"])
        test = dist < 4.
        off = dist > 9.
        row = {"visit":v, "mjd":float(ep.loc[v,"mjd"])}
        for suffix, stack in zip(["all","odd","even"], [stacks[v], *splits[v]]):
            residual = stack-background
            j,i = np.unravel_index(np.nanargmax(np.where(test,residual,np.nan)), xx.shape)
            row["peak_dist_arcsec_" + suffix] = float(dist[j,i])
            row["local_peak_over_offtarget_mad_" + suffix] = float(
                (residual[j,i]-np.nanmedian(residual[off]))/mad(residual[off]))
        panels.append(stacks[v]-background)
        results.append(row)
    pd.DataFrame(results).to_csv(BASE/"image-checks.csv", index=False, float_format="%.5f")
    fig, axes = plt.subplots(len(CHOICES),1, figsize=(5.3,8), sharex=True, sharey=True)
    for ax,v,D in zip(axes,CHOICES,panels):
        ax.imshow(D,origin="lower",extent=[x[0]+SCALE/2,x[-1]-SCALE/2,y[0]-SCALE/2,y[-1]+SCALE/2],
                  cmap="gray",vmin=-3*mad(D),vmax=7*mad(D),aspect="equal")
        ax.add_patch(plt.Circle((ep.loc[v,"x"],ep.loc[v,"y"]),6.5,color="#ff9b30",fill=False,lw=1.3))
        ax.text(0.015,0.85,f'MJD {ep.loc[v,"mjd"]:.0f}',transform=ax.transAxes,color="white",fontsize=9,
                bbox=dict(facecolor="black",alpha=.6,edgecolor="none",pad=2))
        ax.set_yticks([-15,0,15])
    axes[-1].set_xlabel("East offset (arcsec)")
    fig.text(.04,.50,"North offset (arcsec)",ha="center",va="center",rotation=90)
    fig.suptitle("WISE 0855−0714 in NEOWISE W2 difference images",y=.98)
    fig.subplots_adjust(top=.94,bottom=.08,left=.17,right=.98,hspace=.09)
    fig.savefig(BASE/"control-track.png",dpi=130)
    print(pd.DataFrame(results).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
