"""Make same-position LoDeSS photometry and cutout figures in all supplied images.

Usage: python inspect_rows.py DATADIR 9927 2114 8029 8092 13292
"""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.nddata import Cutout2D
from matplotlib.patches import Circle

from photometry import RadioImage

def inspect(root, indices):
    sources = fits.getdata(root/'sources.fits')
    field_map = list(csv.DictReader((root/'source_fields.csv').open()))
    qa = list(csv.DictReader((root/'fields.csv').open()))
    out = []
    meta = []
    for i in indices:
        s = sources[i]
        ra, dec = float(s['RA']), float(s['DEC'])
        sky = SkyCoord(ra*u.deg, dec*u.deg)
        for q in qa:
            sep = sky.separation(SkyCoord(float(q['ra_deg'])*u.deg, float(q['dec_deg'])*u.deg)).deg
            if sep > 10:
                continue
            path = root/'images'/(q['field']+'-I.fits')
            if not path.exists():
                continue
            with RadioImage(path) as radio:
                m = radio.measure(ra, dec)
                if m is None:
                    continue
                stamp = Cutout2D(radio.data, sky, (4*u.arcmin, 4*u.arcmin), wcs=radio.wcs)
                x0,y0 = stamp.to_cutout_position(radio.wcs.world_to_pixel_values(ra, dec))
                hdr = radio.hdus[0].header
                meta.append((stamp.data.copy(), x0,y0, radio.fwhm, m['pixel_noise_jy_beam']))
                out.append(dict(row=i, ra_deg=ra, dec_deg=dec, field=q['field'], quality=q['quality'], centre_deg=sep,
                                parent=field_map[i]['field']==q['field'], date_obs=hdr['DATE-OBS'],
                                catalogue_flux_jy=float(s['Total_flux']), catalogue_peak_jy=float(s['Peak_flux']),
                                catalogue_err_jy=float(s['E_Total_flux']), catalogue_rms_jy=float(s['Isl_rms']), **m))
    if not out:
        return
    with (root/'row_photometry.csv').open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=out[0]);w.writeheader();w.writerows(out)
    n=len(out)
    fig, axes=plt.subplots((n+2)//3,3,figsize=(10.5,3.7*((n+2)//3)),squeeze=False)
    for ax in axes.flat: ax.set_axis_off()
    for ax,r,(image,x,y,fwhm,rms) in zip(axes.flat,out,meta):
        ax.imshow(image/rms,origin='lower',vmin=-3,vmax=8,cmap='RdBu_r',interpolation='nearest')
        ax.add_patch(Circle((x,y),fwhm/2,fill=False,color='lime',lw=1))
        ax.set_title(f"Row {r['row']}, {r['field']}\n{r['psf_jy']*1000:.0f} ± {r['psf_noise_jy']*1000:.0f} mJy",fontsize=10)
    fig.suptitle('LoDeSS 4 arcmin cutouts at catalogue coordinates · colour range −3 to 8 local rms\nGreen circle is the 45 arcsec beam. Errors are empirical local PSF noise only.')
    fig.tight_layout(rect=(0,0,1,.94))
    fig.savefig(root/'row_cutouts.png',dpi=140)
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('rows',nargs='+',type=int)
    args=p.parse_args();inspect(args.root,args.rows)
