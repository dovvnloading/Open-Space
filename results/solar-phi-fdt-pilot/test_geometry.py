"""Synthetic WCS, projection and differential-rotation checks; no images."""
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import HeliographicStonyhurst
from geometry import to_surface,mu,transport,dust_center_zero_based

h=fits.Header()
for k,v in {'DATE-BEG':'2026-08-13T16:25:09','DATE-AVG':'2026-08-13T16:27:22','DATE-END':'2026-08-13T16:29:37','TIMESYS':'UTC','BUNIT':'Gauss','DSUN_OBS':4.75e10,'DSUN_AU':.3175,'HGLN_OBS':-89.48,'HGLT_OBS':-15.75,'RSUN_REF':695700000.,'CTYPE1':'HPLN-TAN','CTYPE2':'HPLT-TAN','CUNIT1':'arcsec','CUNIT2':'arcsec','CRPIX1':960.,'CRPIX2':960.,'CRVAL1':0.,'CRVAL2':0.,'CDELT1':3.6,'CDELT2':3.6,'PXBEG1':65,'PXBEG2':65,'NAXIS1':1920,'NAXIS2':1920,'FILENAME':'synthetic'}.items():h[k]=v
c=to_surface(h,959.,959.);m=mu(h,c); assert abs(m-1)<1e-10
surface=c.transform_to(HeliographicStonyhurst(obstime=h['DATE-AVG']))
assert abs(surface.lon.to_value(u.deg)-h['HGLN_OBS'])<1e-7
assert abs(surface.lat.to_value(u.deg)-h['HGLT_OBS'])<1e-7
x,y,visible=transport(h,h,959.,959.);assert visible and np.hypot(x-959,y-959)<1e-6
c2=to_surface(h,1059.,959.);assert abs(c2.Tx.to_value(u.arcsec)-360.)<.01
assert .95 < mu(h,c2) < 1
x,y,visible=transport(h,h,1059.,959.);assert visible and np.hypot(x-1059,y-959)<1e-6
hrot=h.copy();hrot['PC1_1']=0.;hrot['PC1_2']=-1.;hrot['PC2_1']=1.;hrot['PC2_2']=0.
cp=to_surface(hrot,1059.,959.);assert abs(cp.Ty.to_value(u.arcsec)-360.)<.01
x,y,visible=transport(hrot,hrot,1059.,959.);assert visible and np.hypot(x-1059,y-959)<1e-6
hnext=h.copy();hnext['DATE-AVG']='2026-08-15T16:27:22'
xp,yp,visible=transport(h,hnext,959.,959.);assert visible
xr,yr,visible=transport(hnext,h,xp,yp);assert visible and np.hypot(xr-959,yr-959)<1e-4
h2=h.copy();h2['HGLN_OBS']=h['HGLN_OBS']+180
_,_,visible=transport(h,h2,959.,959.);assert not visible
c3=to_surface(h,10000.,959.);assert not np.isfinite(c3.distance.to_value(u.m))
x,y=dust_center_zero_based(h);assert 1200<x<1300 and 900<y<1000
print('synthetic center/off-axis, rotated WCS, two-day return, far side, off disk and dust indexing pass')
