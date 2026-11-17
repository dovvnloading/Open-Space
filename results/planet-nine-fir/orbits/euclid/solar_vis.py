"""CALSPEC Sun through the Euclid/VIS.vis photon-counting throughput."""
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.io.votable import parse_single_table

HERE=Path(__file__).resolve().parent
FILES={
    'sun_reference_stis_002.fits':'https://ssb.stsci.edu/trds/calspec/sun_reference_stis_002.fits',
    'Euclid_VIS.xml':'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=Euclid/VIS.vis',
}

def main():
    data={}
    receipts={}
    for name,url in FILES.items():
        d=(HERE/'cache'/name)
        d.parent.mkdir(exist_ok=True)
        if not d.exists(): d.write_bytes(urlopen(url,timeout=45).read())
        data[name]=d.read_bytes()
        receipts[name]=dict(url=url,bytes=len(data[name]),sha256=hashlib.sha256(data[name]).hexdigest())
    with fits.open(io.BytesIO(data['sun_reference_stis_002.fits'])) as f:
        sun=f[1].data
        wave=np.array(sun['WAVELENGTH'],dtype=float)
        flux=np.array(sun['FLUX'],dtype=float) # erg cm^-2 s^-1 Angstrom^-1 at 1 au
    band=parse_single_table(io.BytesIO(data['Euclid_VIS.xml'])).array
    w=np.asarray(band['Wavelength']); tr=np.asarray(band['Transmission'])
    # Photon integral for a source and for a constant 3631 Jy AB reference.
    flam=np.interp(w,wave,flux)
    fzero=3631e-23*2.99792458e18/w**2
    mag=-2.5*np.log10(np.trapezoid(flam*tr*w,w)/np.trapezoid(fzero*tr*w,w))
    rows=list(csv.DictReader((HERE.parent/'brightness.csv').open()))
    m=np.array([float(r['r_AB'])+mag+26.93 for r in rows])
    out=dict(sun_VIS_AB_at_1au=float(mag),grey_VIS_AB_model_range=[float(m.min()),float(m.max())],
             inputs=receipts)
    (HERE/'solar_vis.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out)

if __name__=='__main__':main()
