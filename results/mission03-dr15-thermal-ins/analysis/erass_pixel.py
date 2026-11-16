"""Read the eRASS1 pixel upper limit from an official sky tile (no web API required)."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy_healpix import HEALPix
def sample(path,ra,dec):
    ipix=int(HEALPix(nside=65536,order='nested').lonlat_to_healpix(ra*u.deg,dec*u.deg))
    with fits.open(path) as f:
        h=f[1]
        matches=np.flatnonzero(h.data['healpix']==ipix)
        if len(matches)!=1:raise ValueError(f'Expected one pixel; found {len(matches)}')
        row=h.data[matches[0]]
        return {'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'ra':ra,'dec':dec,'ecf':h.header['ECF'],'pixel':{k:row[k].item() for k in h.columns.names}}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('table',type=Path);p.add_argument('ra',type=float);p.add_argument('dec',type=float)
    a=p.parse_args();print(json.dumps(sample(a.table,a.ra,a.dec),indent=2))
