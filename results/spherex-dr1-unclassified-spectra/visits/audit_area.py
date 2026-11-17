"""Check detector coordinates and visit contrast sensitivity to local angular pixel area.

Run after measure.py; takes its output directory and its cached cutout directory.
The finite-difference WCS Jacobian is a local approximation to the solid angle.
It does not supply the mission PSF/aperture correction.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from measure import analyse


def local_area(w, x, y):
    """Sky solid angle in steradians per pixel from opposite mid-edge vectors."""
    rd = np.deg2rad(w.all_pix2world([[x+.5,y],[x-.5,y],[x,y+.5],[x,y-.5]], 0))
    a,b = rd.T
    unit = np.column_stack([np.cos(b)*np.cos(a), np.cos(b)*np.sin(a), np.sin(b)])
    return np.linalg.norm(np.cross(unit[0]-unit[1], unit[2]-unit[3]))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('visits', type=Path)
    ap.add_argument('cache', type=Path)
    args = ap.parse_args()
    reports=[]
    for path in sorted(args.visits.glob('J*.csv')):
        summary = json.loads(path.with_suffix('.json').read_text())
        rows=list(csv.DictReader(path.open()))
        areas, offsets = [], []
        for r in rows:
            r['mjd']=float(r['mjd'])
            if r['status'] != 'ok':
                continue
            hpath = args.cache / (hashlib.sha256(r['url'].encode()).hexdigest()+'.fits')
            h = fits.getheader(hpath,'IMAGE')
            cel = WCS(h)
            x,y = cel.all_world2pix([[summary['ra'],summary['dec']]], 0)[0]
            # Alternate A is explicitly named "0-based active pixel coordinates".
            det = WCS(h,key='A')
            det.sip = None
            dx,dy = det.all_pix2world([[x,y]],0)[0]
            offsets.append(max(abs(dx-float(r['det_x'])),abs(dy-float(r['det_y']))))
            sr = local_area(cel,x,y)
            areas.append(sr * (180*3600/np.pi)**2)
            for name in ('flux','sigma'):
                r[name]=float(r[name])*1e9*sr   # native aperture sum -> approximate mJy
            r['wave_um']=float(r['wave_um'])
        groups=analyse(rows)
        out=dict(target=path.stem, usable=len(areas), min_sq_arcsec=min(areas,default=None),
                 max_sq_arcsec=max(areas,default=None), max_active_pixel_error=max(offsets,default=None))
        for release in ('dr1_all','qr3_all'):
            out[release+'_raw_score']=summary['groups'][release].get('signed_score')
            out[release+'_area_score']=groups[release].get('signed_score')
        reports.append(out)
    if reports:
        with (args.visits/'area_check.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(reports[0]))
            writer.writeheader(); writer.writerows(reports)


if __name__=='__main__':
    main()
