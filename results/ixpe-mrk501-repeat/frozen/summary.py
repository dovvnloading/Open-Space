"""Produce only the authorized baseline summary, suitable for later review.

This script intentionally has no option for selecting a holdout observation.
"""
import json
from astropy.coordinates import SkyCoord
from astropy import units as u
from centroid_rule import baseline
from development import primary,CENTER

if __name__=='__main__':
    pos=baseline();center=SkyCoord(float(pos['ra_deg'])*u.deg,float(pos['dec_deg'])*u.deg)
    print(json.dumps({'centroid':{'ra_deg':float(pos['ra_deg']),'dec_deg':float(pos['dec_deg']),'offset_arcsec':pos['offset_arcsec']},'nominal':primary(CENTER),'xy_centroid':primary(center)},indent=2))
