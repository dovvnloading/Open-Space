"""Convert conditional point-source limits to sizes for specified geometric albedos."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from astropy.io import fits
from astropy.io.votable import parse_single_table

AU_KM=149597870.7
URLS={
    'sun':'https://ssb.stsci.edu/trds/calspec/sun_reference_stis_002.fits',
    'F845M':'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=HST/WFC3_UVIS2.F845M',
    'F162M':'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=JWST/NIRCam.F162M',
    'F140M':'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=JWST/NIRCam.F140M',
}


def calibrate(cache):
    cache.mkdir(exist_ok=True,parents=True);paths={};hashes={}
    for band,url in URLS.items():
        path=cache/(band+'_throughput.xml' if band!='sun' else 'sun_reference_stis_002.fits')
        if not path.exists():
            with urlopen(url,timeout=120) as response:
                path.write_bytes(response.read())
        paths[band]=path;hashes[band]=hashlib.sha256(path.read_bytes()).hexdigest()
    with fits.open(paths['sun']) as h:
        wl,flam=h[1].data['WAVELENGTH'],h[1].data['FLUX']
    result={}
    for name,path in paths.items():
        if name=='sun':continue
        table=parse_single_table(path).to_table()
        wave=np.array(table['Wavelength'],dtype=float); trans=np.array(table['Transmission'],dtype=float)
        if wave[0]<wl.min() or wave[-1]>wl.max():raise ValueError(name+' solar coverage')
        solar=np.interp(wave,wl,flam)
        mean_fnu=np.trapezoid(solar*trans*wave,wave)/(2.99792458e18*np.trapezoid(trans/wave,wave))
        result[name]=float(-2.5*np.log10(mean_fnu)-48.60)
    return dict(solar_AB_at_1au=result,urls=URLS,sha256=hashes)


def diameter(magnitude,sun_mag,distance_au,heliocentric_au,albedo,phase_function=1.):
    return 2*AU_KM*distance_au*heliocentric_au*10**(-.2*(magnitude-sun_mag))/np.sqrt(albedo*phase_function)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path)
    args=ap.parse_args();print(json.dumps(calibrate(args.cache),indent=2))
