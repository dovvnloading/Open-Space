"""Equivalent spherical diameters for reflected-light NIR flux and same-band albedo."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
from astropy.io import fits, votable
import numpy as np

SOLAR = 'https://archive.stsci.edu/hlsps/reference-atlases/cdbs/calspec/sun_reference_stis_002.fits'
FILTER = 'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=JWST/NIRCam.{}'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cube',type=Path)
    p.add_argument('cache',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--band',required=True)
    p.add_argument('--flux-jy',type=float,required=True)
    p.add_argument('--phase-factor',type=float,default=1)
    args = p.parse_args()
    filt = args.band.upper()
    sources = {}
    for name,url in [('solar.fits',SOLAR),(f'{filt.lower()}_filter.xml',FILTER.format(filt))]:
        path = args.cache/name
        if not path.exists():
            args.cache.mkdir(parents=True,exist_ok=True)
            path.write_bytes(urlopen(url,timeout=45).read())
        sources[name]=dict(url=url,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    solar = fits.getdata(args.cache/'solar.fits',1)
    t = votable.parse_single_table(args.cache/f'{filt.lower()}_filter.xml').array
    wav = np.asarray(t['Wavelength'])
    throughput = np.asarray(t['Transmission'])
    flam = np.interp(wav,solar['WAVELENGTH'],solar['FLUX'])
    # Photon-counting mean f_nu: integral(f_lambda*T*lambda) / (c integral(T/lambda)).
    sun_jy = np.trapezoid(flam*throughput*wav,wav)/np.trapezoid(throughput/wav,wav)/2.99792458e18*1e23
    c=np.load(args.cube)
    r,d = np.mean(c['r_au']),np.mean(c['path_delta'])
    albedos=[.04,.08,.16]
    diam = [2*d*149597870.7*r*np.sqrt(args.flux_jy/(sun_jy*alb*args.phase_factor)) for alb in albedos]
    args.output.write_text(json.dumps(dict(sources=sources,r_au=r,observer_distance_au=d,
        filter=filt,flux_Jy=args.flux_jy,ab_mag=-2.5*np.log10(args.flux_jy/3631),
        solar_flux_at_1au_Jy=sun_jy,phase_factor=args.phase_factor,band_geometric_albedos=albedos,
        diameter_km=diam),indent=2)+'\n')


if __name__=='__main__':
    main()
