"""Recover Pan-STARRS epoch astrometry and compare absolute and star-relative tracks."""
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.time import Time
from archival_images import SKY,PM,HOST,HERE

ENDPOINT="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/detection.csv?objID=86863357822955473"
RECEIPT=HERE/'catalogs/ps1_query.json'
DATA=HERE/'catalogs/l_ps_epochs.csv'
PARALLAX=.138228   # parallax in arcseconds
REF_MJD=60977.41317599537

def parallax_shift(mjd):
    ra,dec=np.deg2rad(HOST)
    pos=get_body_barycentric("earth",Time(mjd,format='mjd')).xyz.to_value('au')
    east=np.array([-np.sin(ra),np.cos(ra),0])
    north=np.array([-np.cos(ra)*np.sin(dec),-np.sin(ra)*np.sin(dec),np.cos(dec)])
    return -PARALLAX*np.array([east@pos,north@pos])

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--refresh',action='store_true')
    arg=parser.parse_args()
    if arg.refresh:
        with urlopen(ENDPOINT,timeout=60) as r:
            data=r.read()
            if not data.startswith(b'objID,'):raise ValueError(data[:100])
            DATA.write_bytes(data)
            RECEIPT.write_text(json.dumps(dict(url=ENDPOINT,retrieved_utc=datetime.now(timezone.utc).isoformat(),
                 bytes=len(data),sha256=hashlib.sha256(data).hexdigest()),indent=2)+'\n')
    if hashlib.sha256(DATA.read_bytes()).hexdigest()!=json.loads(RECEIPT.read_text())['sha256']:raise ValueError(DATA)
    detections=list(csv.DictReader(DATA.open()))
    refpar=parallax_shift(REF_MJD)
    out=[]
    for d in sorted(detections,key=lambda d:float(d['obsTime'])):
        t=float(d['obsTime'])
        delta=np.array([(float(d['ra'])-SKY[0])*np.cos(np.deg2rad(SKY[1])),float(d['dec'])-SKY[1]])*3600
        track=np.array(PM)*(t-REF_MJD)/365.25 + parallax_shift(t)-refpar
        out.append(dict(detectID=d['detectID'],mjd=t,date=Time(t,format='mjd').iso[:10],filterID=d['filterID'],
             psfQf=d['psfQf'],ra=d['ra'],dec=d['dec'],fixed_resid_arcsec=np.linalg.norm(delta),
             comoving_resid_arcsec=np.linalg.norm(delta-track),comoving_prediction_shift_arcsec=np.linalg.norm(track),
             track_east_arcsec=track[0],track_north_arcsec=track[1]))
    with (HERE/'ps1_motion.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(out[0]));w.writeheader();w.writerows(out)
    print("detections",len(out),"nights",len(set(v['date'] for v in out)))
    for key in ('fixed_resid_arcsec','comoving_resid_arcsec','comoving_prediction_shift_arcsec'):
        a=[v[key] for v in out];print(key, min(a),np.median(a),max(a))

if __name__=='__main__':main()
