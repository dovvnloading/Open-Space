"""Extract and measure original RACS HiPS samples at the position of V1703 Sco."""
import argparse
import csv
import hashlib
import io
from pathlib import Path
import sys
from urllib.request import urlopen

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parent
NOVA = SkyCoord(267.5802083, -33.6520278, unit="deg")
TILE = "https://casda.csiro.au/hips/RACS/{}/I/Norder8/Dir450000/Npix459508.fits"


def cutout(hdu):
    x, y = WCS(hdu.header).world_to_pixel(NOVA)
    ix, iy = round(float(x)), round(float(y))
    x0, y0 = ix-52, iy-52
    pixels = hdu.data[y0:iy+53,x0:ix+53].copy()
    header = hdu.header.copy()
    header['CRPIX1'] -= x0
    header['CRPIX2'] -= y0
    return fits.PrimaryHDU(data=pixels, header=header)


def measure(path):
    with fits.open(path) as data:
        wcs=WCS(data[0].header)
        image=data[0].data.astype(float)*1000
    x,y=wcs.world_to_pixel(NOVA)
    i,j=np.indices(image.shape)
    radius=np.hypot(j-x,i-y)
    sky=image[(radius>10)&(radius<42)&np.isfinite(image)]
    median=np.median(sky)
    rms=1.482602218*np.median(np.abs(sky-median))
    forced=image[round(float(y)),round(float(x))]
    py,px=np.unravel_index(np.nanargmax(np.where(radius<8,image,np.nan)),image.shape)
    peak=wcs.pixel_to_world(px,py)
    return dict(central_mJy_beam=round(forced,3),background_mJy_beam=round(median,3),
                local_rms_mJy_beam=round(rms,3),snr=round((forced-median)/rms,2),
                local_peak_mJy_beam=round(image[py,px],3),
                local_peak_separation_as=round(NOVA.separation(peak).arcsecond,2),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch',action='store_true',help='Recreate local native-pixel subimages from CSIRO tiles')
    args=parser.parse_args()
    rows=[]
    for band,frequency in [('low',887.5),('mid',1367.5)]:
        url=TILE.format(band)
        path=ROOT/'cutouts'/f'V1703-Sco-RACS-{band}.fits'
        if args.fetch:
            with urlopen(url,timeout=90) as response:
                data=response.read()
            with fits.open(io.BytesIO(data)) as hdus:
                cutout(hdus[0]).writeto(path,overwrite=True)
        rows.append(dict(band=band,frequency_MHz=frequency,original_tile_url=url,**measure(path)))
    writer=csv.DictWriter(sys.stdout,list(rows[0]),lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == '__main__':
    main()
