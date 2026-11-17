"""Screen the entire DESI DR1 SpecDis v2.1 table in chunks. Run from this directory."""
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import astropy.coordinates as coord
import astropy.units as u
import numpy as np
from astropy.io import fits

SOURCE = 'https://data.desi.lbl.gov/public/dr1/vac/dr1/mws-specdis/v2.1/iron-yr1-v2.1.fits'
SHA = '25075043f066c27af6ac25edd381cdff20930ed87f9d337b2d9c2df1403ddf1b'
FRAME = coord.Galactocentric(galcen_distance=8.122*u.kpc, z_sun=20.8*u.pc,
    galcen_v_sun=coord.CartesianDifferential([12.9,245.6,7.78]*u.km/u.s))


def cartesian(ra, dec, dist, pmra, pmdec, rv):
    c = coord.SkyCoord(ra=ra*u.deg, dec=dec*u.deg, distance=dist*u.kpc,
        pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr,
        radial_velocity=rv*u.km/u.s).transform_to(FRAME)
    return c.cartesian.xyz.to_value(u.kpc), c.velocity.d_xyz.to_value(u.km/u.s)


def screen(path):
    with open(path, 'rb') as f:
        assert hashlib.file_digest(f,'sha256').hexdigest() == SHA, 'unexpected input bytes'
    counts = Counter()
    with fits.open(path, memmap=True) as hdus, open('desi_screen.csv','w',newline='') as out:
        data = hdus[1].data
        writer = csv.writer(out)
        writer.writerow(data.names + ['v_spec','v_plx','nearby','fast_spec'])
        for offset in range(0,len(data),80000):
            a = data[offset:offset+80000]
            counts['input_rows'] += len(a)
            with np.errstate(invalid='ignore',divide='ignore'):
                quality = ((a['SOURCE_ID']>0) & (a['RUWE']>0) & (a['RUWE']<1.4)
                    & np.isin(a['APS'],[31,95]) & (a['DIST']>0) & (a['DISTERR']>0)
                    & (a['DISTERR']/a['DIST'] < .25) & (a['BINARY_FLAG']==1)
                    & np.isfinite(a['RA']+a['DEC']+a['PMRA']+a['PMDEC']+a['VRAD']))
            a = a[quality]
            counts['astrometry_specdis_quality_rows'] += len(a)
            p = a['PARALLAX']-a['PARALLAX_ZPC']
            reliable = (a['PARALLAX_ERR']>0) & (p>=5*a['PARALLAX_ERR'])
            counts['corrected_plx_snr5_rows'] += int(reliable.sum())
            def speed(distance,mask):
                return np.linalg.norm(cartesian(a['RA'][mask],a['DEC'][mask],distance,a['PMRA'][mask],a['PMDEC'][mask],a['VRAD'][mask])[1],axis=0)
            vs = speed(a['DIST'],slice(None))
            vp = np.full(len(a),np.nan)
            if reliable.any(): vp[reliable] = speed(1/p[reliable],reliable)
            near = reliable & (vs>=550) & (vp>=550)
            extreme = vs>=900
            counts['specdis_speed_ge550_rows'] += int((vs>=550).sum())
            counts['both_speeds_ge550_rows'] += int(near.sum())
            counts['specdis_speed_ge900_rows'] += int(extreme.sum())
            for i in np.flatnonzero(near|extreme):
                writer.writerow(list(a[i])+[vs[i],vp[i],int(near[i]),int(extreme[i])])
            counts['union_rows'] += int((near|extreme).sum())
    Path('desi_depth.json').write_text(json.dumps(dict(counts),indent=2)+'\n')


if __name__ == '__main__':
    screen(sys.argv[1])
