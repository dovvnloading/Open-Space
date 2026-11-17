"""Verify the extreme-speed subset of the nominal inverse-parallax sensitivity check."""
from pathlib import Path
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord, Galactocentric, CartesianDifferential
import astropy.units as u
d=Path(__file__).resolve().parent
s=pd.read_csv(d/'fast6d.csv',dtype={'targetid':str,'source_id':str})
g=SkyCoord(ra=s.ra.to_numpy()*u.deg,dec=s.dec.to_numpy()*u.deg,distance=(1/s.parallax.to_numpy())*u.kpc,pm_ra_cosdec=s.pmra.to_numpy()*u.mas/u.yr,pm_dec=s.pmdec.to_numpy()*u.mas/u.yr,radial_velocity=s.vrad.to_numpy()*u.km/u.s)
frame=Galactocentric(galcen_distance=8.122*u.kpc,z_sun=20.8*u.pc,galcen_v_sun=CartesianDifferential([11.1,245.6,7.25]*u.km/u.s))
s['v_gc']=g.transform_to(frame).velocity.norm().to_value(u.km/u.s)
s['abs_delta_rv']=abs(s.vrad-299792.458*s.rr_z)
s.loc[s.v_gc>1000,['source_id','targetid','v_gc','vrad','vrad_err','rr_z','abs_delta_rv']].to_csv(d/'parallax_extreme.csv',index=False)
