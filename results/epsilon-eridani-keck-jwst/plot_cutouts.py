"""Plot the two sources in a reference-star frame after running reproduce.py."""
import argparse
from pathlib import Path
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from matplotlib.colors import AsinhNorm
from scipy.ndimage import map_coordinates

from reproduce import EAST_NORTH_SCALE, REFERENCE_SKY, keck_centroid, keck_stack, jw_centroid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    # Use the manifest's five frames, even if the cache contains other Keck products.
    import csv
    with (Path(__file__).parent / "manifest.csv").open() as f:
        manifest = list(csv.DictReader(f))
    keck = [args.cache/r["filename"] for r in manifest if r["instrument"]=="Keck/NIRC2"]
    aligned, scale = keck_stack(keck)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        kimage = np.nanmedian(aligned, axis=0)
    x0, y0 = keck_centroid(kimage, (158,141))
    # Axes are offsets from field source #1. All three cutouts are on the same grid.
    east = np.linspace(-2.25,-4.75,180)
    north = np.linspace(.95,3.45,180)
    ee, nn = np.meshgrid(east,north)
    stamps = [map_coordinates(kimage,[y0-ee/scale,x0-nn/scale],order=1)]
    titles = ["Keck K′  ·  2001 Dec 1"]
    for r in manifest:
        if r["instrument"] != "JWST/NIRCam":
            continue
        with fits.open(args.cache/r["filename"]) as h:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore",FITSFixedWarning)
                wcs = WCS(h["SCI"].header)
            arr = h["SCI"].data
            star, _ = jw_centroid(arr,wcs.all_world2pix(*REFERENCE_SKY,0),9)
            sky = np.array(wcs.all_pix2world(*star,0))
            xx, yy = wcs.all_world2pix(sky[0]+ee/EAST_NORTH_SCALE[0],sky[1]+nn/EAST_NORTH_SCALE[1],0)
            stamps.append(map_coordinates(arr,[yy,xx],order=1))
            titles.append("JWST "+h[0].header["FILTER"]+"  ·  2023 Aug 8")
    fig, axes = plt.subplots(1,3,figsize=(9.3,3.5),sharey=True)
    for ax, arr, title in zip(axes,stamps,titles):
        lo, hi = np.nanpercentile(arr,[4,99.3])
        norm = AsinhNorm(linear_width=.25*(hi-lo),vmin=lo,vmax=hi)
        ax.imshow(arr,origin="lower",extent=[east[0],east[-1],north[0],north[-1]],
                  cmap="gray",norm=norm,interpolation="nearest")
        ax.plot(-3.5,2.2,marker="o",mfc="none",mec="#edaa30",ms=16,mew=1.1)
        ax.set_title(title,fontsize=10)
        ax.set_xticks([-2.5,-3.5,-4.5]);ax.set_yticks([1,2,3]);ax.tick_params(labelsize=8)
    axes[0].set_ylabel("North offset (arcsec)")
    fig.supxlabel("East offset from the reference star (arcsec)",y=.11,fontsize=10)
    fig.text(.5,.015,"Circle marks the published Keck offset. Each panel uses its own intensity scale.",
             ha="center",fontsize=9)
    fig.tight_layout(rect=[0,.14,1,1])
    fig.savefig(args.output,dpi=180,bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
