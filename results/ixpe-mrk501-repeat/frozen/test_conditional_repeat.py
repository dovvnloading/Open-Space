"""Synthetic conditional statistic and GTI-overlap gate; no holdout input."""
import numpy as np
from conditional_repeat import compare,ensure_disjoint,validate_gti
from development import require_region_in_map
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

def pair(q,u):
    x={'q':q,'u':u,'cov_qu':[[.01,0],[0,.01]]}
    return {'nominal':x,'xy_centroid':x}
a=compare(pair(0,0),pair(0,0));assert a['results']['nominal']['conditional_p']==1.
b=compare(pair(0,0),pair(1,1));assert b['results']['nominal']['screening_crossing']
assert abs(b['results']['nominal']['chi2']-100.)<1e-10
ensure_disjoint([(0,1)],[(1,2)])
try:ensure_disjoint([(0,1)],[(.5,2)])
except ValueError:pass
else:raise AssertionError('overlap not rejected')
try:compare(pair(0,0),{'nominal':{'q':float('nan'),'u':0,'cov_qu':[[1,0],[0,1]]},'xy_centroid':pair(0,0)['xy_centroid']})
except ValueError:pass
else:raise AssertionError('NaN accepted')
print('synthetic comparison, GTI and finite gates pass')

validate_gti(np.array([0.,2.]),np.array([1.,3.]))
for starts,stops in [(np.array([0.]),np.array([np.nan])),(np.array([np.nan]),np.array([1.]))]:
    try:validate_gti(starts,stops)
    except ValueError:pass
    else:raise AssertionError('non-finite GTI accepted')
w=WCS(naxis=2);w.wcs.crpix=[300,300];w.wcs.cdelt=[-.00072222,.00072222];w.wcs.crval=[253.46758,39.76017];w.wcs.ctype=['RA---TAN','DEC--TAN']
center=SkyCoord(253.46758*u.deg,39.76017*u.deg)
require_region_in_map(w,(600,600),center)
edge=SkyCoord(*w.all_pix2world(530,300,0),unit='deg')
try:require_region_in_map(w,(600,600),edge)
except ValueError:pass
else:raise AssertionError('clipped background annulus accepted')
print('non-finite GTI and clipped-annulus gates pass')
