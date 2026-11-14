from pathlib import Path
import json
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation,SkyCoord
from astropy import units as u
from astropy.utils import iers
iers.conf.auto_download=False
D=Path(__file__).resolve().parent
SKY=SkyCoord(21.06275*u.deg,-73.32594444*u.deg)
PARAMS=json.loads((D/'inputs/ogle-model.json').read_text())['params']
PERIOD=1.1696719
def hjd(mjd):
    # Barycentric/heliocentric difference of orbit around Earth is small at this cadence.
    t=Time(mjd,format='mjd',scale='utc',location=EarthLocation.from_geocentric(0.,0.,0.,u.m))
    return t.jd+t.light_travel_time(SKY,kind='heliocentric').to_value(u.day)
def shape(jd):
    p=PARAMS
    year=(jd-2456000.)/365.25
    phase=((jd-2460000.0369)/PERIOD)%1
    shift=p[10]*np.cos(2*np.pi*year/p[12])+p[11]*np.sin(2*np.pi*year/p[12])
    c1=p[2]+p[8]*year-shift
    c2=p[3]+p[8]*year+shift
    return phase,((phase-c1+.5)%1)-.5,((phase-c2+.5)%1)-.5
