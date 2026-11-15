"""Plot all four epochs on one east/north tangent plane."""
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from scipy.ndimage import map_coordinates
from replay import Image, predicted, ROOT

def main():
    images = [Image(p) for p in sorted((ROOT/"inputs").glob("VLASS*.fits"))]
    ra0, dec0 = predicted(images[0].date)
    grid = np.linspace(-8,8,161)
    e, n = np.meshgrid(grid,grid)
    ra, dec = ra0+e/3600/np.cos(np.deg2rad(dec0)), dec0+n/3600
    fig, axes = plt.subplots(1,4,figsize=(10.8,3.4),layout="constrained",sharey=True)
    for im, ax in zip(images,axes):
        x, y = im.wcs.world_to_pixel_values(ra,dec)
        data = map_coordinates(im.image,[y,x],order=1,mode="nearest")
        raster=ax.imshow(data,origin="lower",extent=[-8,8,-8,8],cmap="magma",vmin=-.25,vmax=.7)
        p=predicted(im.date)
        ep=(p[0]-ra0)*3600*np.cos(np.deg2rad(dec0))
        npix=(p[1]-dec0)*3600
        ax.scatter([0],[0],s=80,marker="+",color="white",lw=1.2)
        ax.scatter([ep],[npix],s=115,marker="o",facecolors="none",edgecolors="#55e6e0",lw=1.35)
        ax.add_patch(Ellipse((5.6,-5.6),im.axes[0],im.axes[1],angle=90-np.rad2deg(im.pa),fc=".8",ec="none"))
        ax.set_title(im.date[:10],fontsize=11)
        ax.set_xlim(8,-8)
        ax.set_xticks([6,3,0,-3,-6]);ax.set_yticks([-6,-3,0,3,6])
    axes[0].set_ylabel("North offset (arcsec)")
    fig.supxlabel("East offset from 2017 infrared prediction (arcsec)",fontsize=10)
    fig.suptitle("WISE J0458+6434AB in the four VLASS visits",x=.01,ha="left",fontsize=14)
    fig.colorbar(raster,ax=axes,label="mJy per beam",shrink=.8)
    fig.savefig(ROOT/"four-epochs.png",dpi=165)
if __name__=="__main__":
    main()
