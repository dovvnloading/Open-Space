"""Aperture photometry of the 10 original PS1 i-band warps covering both contrasts."""
import argparse
import csv
import hashlib
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

# WD and two nearest persistent neighbors in the MAST DR2 mean catalog.
STARS = [
    ("wd", 140.96429, 73.43997),
    ("neighbor86", 140.95814, 73.44155),
    ("neighbor96", 140.95833, 73.43786),
]
ARCHIVE="https://ps1images.stsci.edu/cgi-bin/"
RA, DEC = 140.9643333333, 73.4399111111


def measure(im, x, y, radius):
    yy, xx = np.indices(im.shape)
    r = np.hypot(xx-x, yy-y)
    # 0.25 arcsec/pixel. These local annuli lie inside the WD–neighbor spacing.
    annulus = im[(r>=17) & (r<=24)]
    sky = np.nanmedian(annulus)
    mask = r < radius
    pix = im[mask] - sky
    if not np.isfinite(pix).all() or np.isfinite(annulus).mean() < .9:
        return np.nan, np.nan
    flux = pix.sum()
    offset = np.hypot(np.sum(pix*(xx[mask]-x)), np.sum(pix*(yy[mask]-y))) / flux
    return float(flux), float(offset*.25)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a=p.parse_args()
    a.cache.mkdir(exist_ok=True,parents=True)
    a.out.mkdir(exist_ok=True,parents=True)
    query = ARCHIVE + "ps1filenames.py?" + urlencode(dict(ra=RA,dec=DEC,type="warp",filters="iz",format="csv"))
    manifest = a.cache/"warp.csv"
    if not manifest.exists():
        manifest.write_bytes(urlopen(query,timeout=150).read())
    lines = manifest.read_text().splitlines()
    records = []
    for line in lines[1:]:
        d=dict(zip(lines[0].split(),line.split()))
        mjd=float(d["mjd"])
        if d["filter"]!="i" or int(mjd) not in (56289,56319,56676,56729) or int(d["badflag"]):
            continue
        url=ARCHIVE+"fitscut.cgi?"+urlencode(dict(ra=RA,dec=DEC,size=240,format="fits",red=d["filename"]))
        file=a.cache/"warps"/d["shortname"]
        file.parent.mkdir(exist_ok=True,parents=True)
        if not file.exists():
            file.write_bytes(urlopen(url,timeout=150).read())
        im, header = fits.getdata(file,header=True)
        coords=WCS(header,naxis=2).all_world2pix([(ra,dec) for _,ra,dec in STARS],0)
        for radius in [6,8,10,12]:
            result = dict(mjd=mjd, aperture_arcsec=radius*.25, filename=d["shortname"], sha256=hashlib.sha256(file.read_bytes()).hexdigest())
            for (name,_,_), (x,y) in zip(STARS,coords):
                result[name], result[name+"_centroid_offset_arcsec"]=measure(im,x,y,radius)
            for name,_,_ in STARS[1:]:
                result["mag_rel_"+name]=-2.5*np.log10(result["wd"]/result[name])
            records.append(result)
    with (a.out/"pixel-photometry.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)


if __name__ == "__main__":
    main()
