"""Re-measure U in the public, calibrated ALMA target continuum maps."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import tarfile
from urllib.request import urlopen

import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent


def measure(item, cache):
    path = cache / ('almaextra-' + item['short_id'] + '.tar')
    if not path.exists():
        with urlopen(item['url'], timeout=60) as response:
            path.write_bytes(response.read())
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == item['sha256']
    with tarfile.open(fileobj=io.BytesIO(data)) as archive:
        raw = archive.extractfile(item['member']).read()
    assert hashlib.sha256(raw).hexdigest() == item['member_sha256']
    hdu = fits.open(io.BytesIO(raw))[0]
    h = hdu.header
    assert h['BUNIT'].lower() == 'jy/beam'
    image = 1000*hdu.data.squeeze().astype(float)
    w = WCS(h).celestial
    ra = 15*(14+39/60+28.505/3600)
    dec = -(60+49/60+52.50/3600)
    xp, yp = w.world_to_pixel_values(ra,dec)
    yy,xx = np.indices(image.shape)
    east = (xx-xp)*h['CDELT1']*3600
    north = (yy-yp)*h['CDELT2']*3600
    radius = np.hypot(east,north)
    fwhm, minor = h['BMAJ']*3600, h['BMIN']*3600
    pa = np.radians(h['BPA'])
    near = (radius < 2.1*fwhm) & np.isfinite(image)
    es, ns = east[near],north[near]
    def residual(parameters):
        flux,x,y,sky,gx,gy=parameters
        dx,dy=es-x,ns-y
        u,v=dx*np.sin(pa)+dy*np.cos(pa), dx*np.cos(pa)-dy*np.sin(pa)
        model=np.exp(-4*np.log(2)*((u/fwhm)**2+(v/minor)**2))
        return sky+gx*es+gy*ns+flux*model-image[near]
    p=least_squares(residual,[max(np.max(image[near]),0.01),0,0,0,0,0],
                    bounds=([0,-1.5,-1.5,-np.inf,-np.inf,-np.inf],
                            [np.inf,1.5,1.5,np.inf,np.inf,np.inf])).x
    fixed_beam = residual([1,0,0,0,0,0])-residual([0,0,0,0,0,0])
    fixed_flux = np.linalg.lstsq(np.column_stack([fixed_beam, np.ones_like(es),es,ns]),image[near],rcond=None)[0][0]
    center_x,center_y=h['CRPIX1']-1,h['CRPIX2']-1
    dist=np.hypot((xx-center_x)*h['CDELT1'],(yy-center_y)*h['CDELT2'])*3600
    source_dist=np.hypot((xp-center_x)*h['CDELT1'],(yp-center_y)*h['CDELT2'])*3600
    control=(abs(dist-source_dist)<0.7*fwhm) & (radius>3*fwhm) & np.isfinite(image)
    v=image[control]
    noise=1.4826*np.median(np.abs(v-np.median(v)))
    ra_fit,dec_fit=w.pixel_to_world_values(xp+p[1]/(3600*h['CDELT1']),yp+p[2]/(3600*h['CDELT2']))
    return dict(product=item['short_id'],utc_header=h['DATE-OBS'],ra_deg=float(ra_fit),dec_deg=float(dec_fit),
                fixed_peak_mjy_per_beam=fixed_flux,fixed_peak_to_rms=fixed_flux/noise,offset_east_arcsec=p[1],offset_north_arcsec=p[2],mjy_per_beam=p[0],
                annulus_mad_rms_mjy_per_beam=noise,peak_to_rms=p[0]/noise,
                fwhm_major_arcsec=fwhm,fwhm_minor_arcsec=minor,control_pixels=len(v))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('cache',type=Path)
    args=parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    products=json.loads((HERE/'receipts.json').read_text())['products']
    rows=[measure(i,args.cache) for i in products]
    with (HERE/'alma_pixels.csv').open('w') as f:
        writer=csv.DictWriter(f,list(rows[0])); writer.writeheader(); writer.writerows(rows)


if __name__=='__main__':
    main()
