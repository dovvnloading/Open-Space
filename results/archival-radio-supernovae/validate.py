"""Restoring-beam fits, sky-registered comparison positions and free radio centroids."""
import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from scipy.optimize import least_squares
from scipy.signal import fftconvolve

def fit_field(path, ra, dec):
    with fits.open(path) as hdus:
        header = hdus[0].header
        image = hdus[0].data.squeeze().astype(float) * 1000
    wcs = WCS(header).celestial
    x, y = map(float, wcs.world_to_pixel_values(ra, dec))
    ix, iy = round(x), round(y)
    dy, dx = np.mgrid[-8:9, -8:9]
    dx, dy = dx + ix - x, dy + iy - y
    pa = np.radians(header['BPA'])
    major, minor = header['BMAJ'] * 3600, header['BMIN'] * 3600
    angular = wcs.pixel_scale_matrix * 3600
    def beam(dx, dy):
        east = angular[0, 0] * dx + angular[0, 1] * dy
        north = angular[1, 0] * dx + angular[1, 1] * dy
        u = (east * np.sin(pa) + north * np.cos(pa)) / major
        v = (east * np.cos(pa) - north * np.sin(pa)) / minor
        return np.exp(-4 * np.log(2) * (u*u + v*v))
    basis = np.array([beam(dx, dy), np.ones(dx.shape), dx, dy]).reshape(4, -1).T
    kernel = np.linalg.pinv(basis)[0].reshape(dx.shape)
    filt = fftconvolve(np.nan_to_num(image), kernel[::-1, ::-1], mode='same')
    valid = fftconvolve(np.isfinite(image).astype(float), np.ones(dx.shape), mode='same') > dx.size - .5
    filt[~valid] = np.nan
    target = filt[iy, ix]
    patch = image[iy-8:iy+9, ix-8:ix+9]
    def residual(params):
        amp, u, v, bg, gx, gy = params
        return (amp * beam(dx-u, dy-v) + bg + gx*dx + gy*dy - patch).ravel()
    fit = least_squares(residual, [max(target,.05),0,0,np.nanmedian(patch),0,0],
                        bounds=([0,-3,-3,-np.inf,-np.inf,-np.inf],
                                [np.inf,3,3,np.inf,np.inf,np.inf])).x
    centroid = wcs.pixel_to_world_values(x+fit[1], y+fit[2])
    # Each position lies at the same tangent-plane offset from the optical sky position in all visits.
    dn,de=np.mgrid[-180:181:6,-180:181:6]
    use=np.hypot(de,dn)>=12
    de,dn=de[use],dn[use]
    cx,cy=wcs.world_to_pixel_values(ra+de/(3600*np.cos(np.deg2rad(dec))),dec+dn/3600)
    controls=map_coordinates(filt,[cy-(y-iy),cx-(x-ix)],order=1,cval=np.nan)
    local = controls[np.isfinite(controls)]
    med=np.median(local); noise=1.4826*np.median(abs(local-med))
    return dict(forced=target,scatter=noise,median=med,free=fit[0],
                sep=np.linalg.norm(angular@fit[1:3]),cen_ra=float(centroid[0]),cen_dec=float(centroid[1]),
                beam_major=major,beam_minor=minor,beam_pa=header['BPA'],date=header['DATE-OBS'],
                n_control=len(local),controls_ge=int(sum(local>=target)),
                field_sha256=hashlib.sha256(path.read_bytes()).hexdigest()),controls

def analyze(rows, directory, output):
    measurements, pairs=[],[]
    for name,group in rows.groupby('name'):
        values=[]
        for row in group.sort_values('epoch').itertuples():
            path=directory/f'{name}-{row.epoch}.fits'
            if not path.exists():continue
            m,ctrl=fit_field(path,row.ra,row.dec)
            m.update(name=name,epoch=row.epoch,days=row.mjd-row.peakmjd)
            measurements.append(m);values.append((m,ctrl))
        for a,(ma,ca) in enumerate(values):
            for mb,cb in values[a+1:]:
                diff=cb-ca; diff=diff[np.isfinite(diff)]
                if not len(diff):continue
                change=mb['forced']-ma['forced'];center=np.median(diff)
                scale=1.4826*np.median(abs(diff-center))
                pairs.append(dict(name=name,first=ma['epoch'],second=mb['epoch'],change=change,
                    control_median=center,control_mad=scale,field_z=(change-center)/scale,
                    n=len(diff),abs_changes_ge=int(sum(abs(diff-center)>=abs(change-center)))))
    pd.DataFrame(measurements).to_csv(output/'beam-fits.csv',index=False,float_format='%.8g')
    pd.DataFrame(pairs).to_csv(output/'field-differences.csv',index=False,float_format='%.8g')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ledger',type=Path);parser.add_argument('images',type=Path)
    parser.add_argument('output',type=Path);args=parser.parse_args()
    analyze(pd.read_csv(args.ledger),args.images,args.output)
