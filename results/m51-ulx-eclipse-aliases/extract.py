"""Extract ACIS photons and chip-specific good time intervals around the two M51 ULXs."""
from pathlib import Path
import argparse
import csv
import gzip
import warnings

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

# Urquhart & Soria (2016), Sec. 1, decimal degrees, ICRS.
TARGETS = [(15*(13+29/60+39.94/3600), 47+12/60+36.6/3600),
           (15*(13+29/60+39.44/3600), 47+12/60+43.3/3600)]
RA, DEC = np.deg2rad(TARGETS[0])
UNIT = np.array([np.cos(RA)*np.cos(DEC), np.sin(RA)*np.cos(DEC), np.sin(DEC)])

def extract(path, out):
    with fits.open(path, memmap=False) as hdus:
        events = hdus['EVENTS']
        h, data = events.header, events.data
        cols = [events.columns.names.index(n)+1 for n in ['x','y']]
        sky = WCS(naxis=2)
        for key, attr in [('TCRVL','crval'),('TCRPX','crpix'),('TCDLT','cdelt'),('TCTYP','ctype')]:
            setattr(sky.wcs, attr, [h[f'{key}{i}'] for i in cols])
        centers = sky.all_world2pix(TARGETS, 1)
        offsets = [(data['x']-x, data['y']-y) for x,y in centers]
        scale = abs(sky.wcs.cdelt[0])*3600
        r1, r2 = [scale*np.hypot(dx,dy) for dx,dy in offsets]
        energy = data['energy']
        keep = ((r1<=7)|(r2<=2.5)) & (energy>=300)&(energy<=8000)
        if not np.any(keep):
            print(h['OBS_ID'], 'no events within target region; nearest event (arcsec)', float(r1.min()))
            return
        chip = int(np.bincount(data['ccd_id'][keep]).argmax())
        if np.any(data['ccd_id'][keep]!=chip):
            raise ValueError(f'multiple CCDs at sources in {path}')
        gti = next(hdu.data for hdu in hdus if hdu.name=='GTI' and hdu.header['CCD_ID']==chip)
        # Chandra orbit adds at most ~0.5 s. Interpolate the geocenter
        # light-time correction over the visit rather than per CCD frame.
        t0, t1 = gti['START'].min(), gti['STOP'].max()
        knots = np.linspace(t0,t1, max(2,int((t1-t0)/3600)+2))
        times = Time(h['MJDREF']+(knots+h.get('TIMEZERO',0))/86400,format='mjd',scale='tt')
        earth = get_body_barycentric('earth', times.tdb).xyz.to_value('au')
        bjd = times.tdb.mjd + (UNIT @ earth)*499.004783836/86400
        def to_mjd(t):
            return np.interp(t,knots,bjd)
        bounds = np.column_stack((gti['START'],gti['STOP']))
        inside = np.zeros(keep.sum(),dtype=bool)
        for start, stop in bounds:
            inside |= (data['time'][keep]>=start)&(data['time'][keep]<=stop)
        inds=np.flatnonzero(keep)[inside]
        obs=h['OBS_ID']
        with gzip.open(out/f'{obs}.csv.gz','wt',newline='') as f:
            w=csv.writer(f); w.writerow(['mjd_tdb_geo','energy_ev','dx_arcsec','dy_arcsec','r1_arcsec','r2_arcsec'])
            for t,e,dx,dy,a,b in zip(to_mjd(data['time'][inds]),energy[inds],offsets[0][0][inds]*scale,offsets[0][1][inds]*scale,r1[inds],r2[inds]):
                w.writerow([f'{t:.9f}',f'{e:.1f}',f'{dx:.4f}',f'{dy:.4f}',f'{a:.4f}',f'{b:.4f}'])
        with (out/f'{obs}-gti.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['obsid','ccd','start_mjd_tdb_geo','stop_mjd_tdb_geo'])
            for lo,hi in to_mjd(bounds):w.writerow([obs,chip,f'{lo:.9f}',f'{hi:.9f}'])

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache',type=Path)
    parser.add_argument('out',type=Path)
    a=parser.parse_args(); a.out.mkdir(exist_ok=True,parents=True)
    for path in sorted(a.cache.glob('*.fits.gz')):
        extract(path,a.out)
