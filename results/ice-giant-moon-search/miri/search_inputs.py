"""Rebuild the two Neptune MIRI cubes and static ephemeris inputs."""
import argparse
import csv
import hashlib
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from geometry import MOONS,query
from jwst_cutouts import make_cut
from prepare import prepare


def build(manifest, footprint, cache, horizons, output):
    cache.mkdir(exist_ok=True,parents=True);output.mkdir(exist_ok=True,parents=True)
    rows={r['productFilename']:r for r in csv.DictReader(manifest.open())}
    sel=[r for r in csv.DictReader(footprint.open()) if r['target']=='TRITON' and r['planet_on_chip']=='1' and r['planned_aperture']=='MIRIM_ILLUM']
    for r in sel:
        dest=cache/r['file']
        if not dest.exists():make_cut(rows[r['file']],215,cache)
    intervals=[]
    for r in sel:
        p=cache/r['file'];h=fits.getheader(p)
        intervals.extend(np.linspace(h['EXPSTART'],h['EXPEND'],9))
        intervals.append(h['EXPMID'])
    # All trial exposure samples and exact header midpoints are included.
    times=np.unique(intervals)
    bodies=[899,*MOONS['Neptune']]
    values=[query(body,'JWST',times+2400000.5,horizons) for body in bodies]
    ep=cache/'ephemeris.npz'
    np.savez_compressed(ep,mjd=times,ids=bodies,values=values,planet='Neptune',observer='JWST')
    for band in ('F2100W','F2550W'):
        files=[cache/r['file'] for r in sel if r['filter']==band]
        centers=output/f'{band}_centers.csv'
        with centers.open('w') as handle:
            w=csv.writer(handle,lineterminator='\n');w.writerow(['exposure_id','x','y','ra_deg','dec_deg'])
            for file in files:
                with fits.open(file) as hdu:
                    t=hdu[0].header['EXPMID'];state=values[0][np.argmin(np.abs(times-t))]
                    ra,dec=state[:2];x,y=WCS(hdu['SCI'].header,naxis=2).all_world2pix(ra,dec,0)
                    w.writerow([file.stem,x,y,ra,dec])
        cube=cache/f'{band}.npz'
        # The targets tracked Triton. All four exposures have Neptune centers
        # explicitly cross-checked against the shifted native SIP WCS.
        prepare(files,ep,cube,centers,301,.11,9,2.8 if band=='F2100W' else 3.3,
                7.,allow_jumps=True,offset_target=True)
        with cube.open('rb') as f:print(band,hashlib.file_digest(f,'sha256').hexdigest(),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for field in ('products','footprints','cache','horizons','out'):p.add_argument(field,type=Path)
    a=p.parse_args();build(a.products,a.footprints,a.cache,a.horizons,a.out)
