"""Compute aperture centroids in celestial coordinates for the 3.4-um galaxy scan."""
import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

RA,DEC=209.033074,30.0813453
ROOT=Path(__file__).resolve().parent

def run(cache):
    records=[]
    with (ROOT/'profile_J135607.csv').open() as f:
        for r in csv.DictReader(f):
            if r['status']!='ok' or not 3.07<=float(r['wave_um'])<=3.73:
                continue
            path=cache/'cutouts'/(hashlib.sha256(r['url'].encode()).hexdigest()+'.fits')
            with fits.open(path) as hd:
                h=hd['IMAGE'].header
                w=WCS(h)
                x,y=w.all_world2pix([[RA,DEC]],0)[0]
                img=hd['IMAGE'].data-hd['ZODI'].data
                yy,xx=np.indices(img.shape)
                radius=np.hypot(xx-x,yy-y)
                mask=(radius<=2.25)
                clean=np.isfinite(img)&((hd['FLAGS'].data&~(1<<21))==0)
                sky=np.median(img[(radius>=4)&(radius<=7)&clean])
                weights=img[mask]-sky
                if np.sum(weights)<=0:continue
                cx,cy=np.sum(xx[mask]*weights)/np.sum(weights),np.sum(yy[mask]*weights)/np.sum(weights)
                cra,cdec=w.all_pix2world([[cx,cy]],0)[0]
                dx=(cra-RA)*np.cos(np.deg2rad(DEC))*3600;dy=(cdec-DEC)*3600
                records.append(dict(release=r['release'],obs_id=r['obs_id'],wave_um=r['wave_um'],
                                    ra_offset_arcsec=round(dx,3),dec_offset_arcsec=round(dy,3),
                                    radial_offset_arcsec=round(np.hypot(dx,dy),3)))
    with (ROOT/'profile_positions.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);run(p.parse_args().cache)
