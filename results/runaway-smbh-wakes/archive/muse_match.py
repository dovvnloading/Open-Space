"""Measure along- and cross-track position for independent MUSE HUDF DR2 spectroscopic redshifts."""
import argparse,hashlib,io,json,tarfile
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
URL='https://amused.univ-lyon1.fr/data/UDF/HUDF/download/UDF_HUDF_Catalogs_FULL.tar.gz'
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('images',type=Path);p.add_argument('results',type=Path)
    p.add_argument('output',type=Path);a=p.parse_args()
    archive=a.images/'muse_catalog.tar.gz'
    if not archive.exists():
        with urlopen(URL,timeout=150) as r:archive.write_bytes(r.read())
    with tarfile.open(archive) as t: payload=t.extractfile('dr2_main_09.fits').read()
    d=fits.open(io.BytesIO(payload))[1].data
    w=WCS(fits.getheader(next(a.images.glob('*f606w_v1_sci.fits')))).celestial
    row=json.loads((a.results/'xdf_w_606_profile.json').read_text())['detections'][59]
    xy=np.array([row['x0'],row['y0']]);v=np.array([row['x1'],row['y1']])-xy;v/=np.linalg.norm(v)
    sky=SkyCoord(d['RA'],d['DEC'],unit='deg')
    sep=sky.separation(SkyCoord(*w.all_pix2world([[row['x'],row['y']]],0)[0],unit='deg')).arcsec
    candidates=[]
    for i in np.argsort(sep):
        if sep[i]>10.8:break
        delta=(np.array(w.world_to_pixel(sky[i]))-xy)*.06
        candidates.append(dict(id=int(d['ID'][i]),raf_id=int(d['RAF_ID'][i]),ra=float(d['RA'][i]),dec=float(d['DEC'][i]),
             z=float(d['Z'][i]),zconf=int(d['ZCONF'][i]),separation_arcsec=sep[i],
             along_arcsec=float(delta@v),cross_arcsec=float(delta@[-v[1],v[0]])))
    result=dict(url=URL,archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),member='dr2_main_09.fits',
                member_sha256=hashlib.sha256(payload).hexdigest(),total_rows=len(d),cone_radius_arcsec=10.8,matches=candidates)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
