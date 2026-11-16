"""Compare aperture flux at native and matched resolution, using official LoTSS DR3 cutouts.

Usage: python resolution_check.py DATA_DIRECTORY COMPARISON_RESULTS_DIRECTORY ROW [ROW ...]
"""
import argparse
from pathlib import Path
from urllib.parse import urlencode

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from matplotlib.patches import Circle
from scipy.ndimage import gaussian_filter

from acquire_lodess import retrieve, digest
from photometry import mad

def aperture(a, x0, y0, pix_arcsec, fwhm, ra=60.):
    yy, xx = np.ogrid[:a.shape[0], :a.shape[1]]
    r = pix_arcsec*np.hypot(xx-x0, yy-y0)
    sky = (r>=90)&(r<=150)
    bg = np.nanmedian(a[sky])
    e = np.exp(-4*np.log(2)*(ra/fwhm)**2)
    return np.nansum(a[r<=ra]-bg)/(np.pi/(4*np.log(2))*(fwhm/pix_arcsec)**2*(1-e)), bg, mad(a[sky])

def compare(root, matches, ids):
    l = pd.read_csv(matches/'match_lotss.csv.gz').set_index('row')
    q = pd.read_csv(root/'source_fields.csv').set_index('row')
    fig, axes = plt.subplots(len(ids),3,figsize=(10,3.6*len(ids)),squeeze=False)
    products = []
    results = []
    (root/'lotss_cutouts').mkdir(exist_ok=True)
    for i, row in enumerate(ids):
        ra, dec = l.loc[row,['lotss_ra','lotss_dec']]
        sky=SkyCoord(ra*u.deg, dec*u.deg)
        pos=sky.to_string(style='hmsdms',sep=':',precision=3)
        for j, resolution in enumerate([45,20,6]):
            if resolution==45:
                f=root/'images'/(q.loc[row,'field']+'-I.fits')
            else:
                suffix='-low' if resolution==20 else ''
                url=f'https://lofar-surveys.org/dr3{suffix}-cutout.fits?'+urlencode({'pos':pos,'size':12})
                f=root/'lotss_cutouts'/f'{row}_{resolution}.fits'
                # Re-download for exact coordinate and URL provenance.
                retrieve(url,f)
                products.append(dict(filename=f.name,sha256=digest(f),url=url))
            with fits.open(f,memmap=True) as h:
                hdr=h[0].header
                w=WCS(hdr).celestial
                stamp=Cutout2D(np.squeeze(h[0].data),sky,(7*u.arcmin,7*u.arcmin),wcs=w)
                arr=stamp.data.astype(float)
                x,y=stamp.wcs.world_to_pixel(sky)
                scale=abs(hdr['CDELT1'])*3600
                beam=hdr['BMAJ']*3600
                f_native,bg,noise=aperture(arr,x,y,scale,beam)
                f45=f_native
                if resolution!=45:
                    sigma=np.sqrt(45**2-beam**2)/(np.sqrt(8*np.log(2))*scale)
                    smoothed=gaussian_filter(arr-bg,sigma=sigma,mode='nearest')*(45/beam)**2
                    f45,_,_=aperture(smoothed,x,y,scale,45)
                results.append(dict(row=row, survey='LoDeSS' if j==0 else 'LoTSS',native_beam_arcsec=beam,
                                    aperture_radius_arcsec=60, native_jy=f_native,matched_45_jy=f45,
                                    native_rms_jy_beam=noise))
                ax=axes[i,j]
                crop=Cutout2D(arr-bg,(x,y),(100*u.arcsec,100*u.arcsec),wcs=stamp.wcs)
                top=max(5*noise,np.nanmax(crop.data)*.65)
                ax.imshow(crop.data,origin='lower',cmap='inferno',interpolation='nearest',vmin=-2*noise,vmax=top)
                ax.add_patch(Circle((.17*crop.data.shape[1],.17*crop.data.shape[0]),beam/(2*scale),fill=False,color='w',lw=1))
                ax.set_title(('LoDeSS 23 MHz' if j==0 else 'LoTSS 144 MHz')+f', {beam:.0f}″ beam\n{f_native:.3f} Jy in 60″ radius',fontsize=10)
                ax.set_xticks([]);ax.set_yticks([])
                if j==0:ax.set_ylabel(f"row {row}",weight='bold')
    fig.suptitle('Common sky footprint (100″ per panel) · individual brightness scales')
    fig.tight_layout(rect=(0,0,1,.94))
    fig.savefig(root/'resolution_check.png',dpi=150);plt.close(fig)
    pd.DataFrame(results).to_csv(root/'resolution_photometry.csv',index=False)
    pd.DataFrame(products).to_csv(root/'resolution_manifest.csv',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('matches',type=Path);p.add_argument('rows',nargs='+',type=int)
    args=p.parse_args();compare(args.root,args.matches,args.rows)
