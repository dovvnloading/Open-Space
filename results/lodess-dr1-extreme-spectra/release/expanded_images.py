"""Reproduce the extended-source environment and aperture sweep for J115948.8+582020."""
import argparse
from pathlib import Path
from urllib.parse import urlencode

import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.patches import Circle
from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pandas as pd

from acquire_lodess import retrieve, digest
from resolution_check import aperture

def inspect(root, out):
    out.mkdir(parents=True,exist_ok=True)
    pos='11:59:48.8 +58:20:20'
    sky=SkyCoord(179.953484*u.deg,58.338758*u.deg)
    urls={'l54':'https://lofar-surveys.org/public/lolss/mosaic-i/p179+57-mosaic.fits'}
    for resolution in [6,20]:
        suffix='-low' if resolution==20 else ''
        urls[f'l144-{resolution}']='https://lofar-surveys.org/dr3'+suffix+'-cutout.fits?'+urlencode({'pos':pos,'size':20})
    records=[]
    for key,url in urls.items():
        f=out/(key+'.fits')
        retrieve(url,f)
        records.append(dict(file=f.name,sha256=digest(f),url=url))
    pd.DataFrame(records).to_csv(out/'j115948_manifest.csv',index=False)
    fig,axes=plt.subplots(1,4,figsize=(14,3.8))
    results=[]
    inputs=[('LoDeSS 23 MHz',root/'images/P181+55-I.fits'),('LoLSS 54 MHz',out/'l54.fits'),
            ('LoTSS 144 MHz',out/'l144-20.fits'),('LoTSS 144 MHz',out/'l144-6.fits')]
    for ax,(lab,f) in zip(axes,inputs):
        with fits.open(f) as h:
            hdr=h[0].header
            a=Cutout2D(h[0].data.squeeze(),sky,(12*u.arcmin,12*u.arcmin),wcs=WCS(hdr).celestial)
            scale=abs(hdr['CDELT1'])*3600
            beam=hdr['BMAJ']*3600
            x,y=a.wcs.world_to_pixel(sky)
            for radius in [30,45,60,90]:
                flux,bg,noise=aperture(a.data,x,y,scale,beam,ra=radius)
                results.append(dict(survey=lab,beam=beam,rad=radius,flux_jy=flux,bg=bg,pixel_rms=noise))
            cut=Cutout2D(a.data-bg,sky,(8*u.arcmin,8*u.arcmin),wcs=a.wcs)
            vx,vy=cut.wcs.world_to_pixel(sky)
            norm=SymLogNorm(linthresh=noise*2,vmin=-noise*3,vmax=np.nanmax(cut.data))
            ax.imshow(cut.data,origin='lower',cmap='inferno',norm=norm,interpolation='nearest')
            ax.plot(vx,vy,marker='+',color='cyan',ms=10,mew=.8)
            ax.add_patch(Circle((vx,vy),22.5/scale,fill=False,color='cyan',lw=.5))
            ax.set_title(f'{lab} · {beam:.0f}″')
            ax.set_xticks([])
            ax.set_yticks([])
    fig.suptitle('J115948.8+582020 and surrounding radio emission · 8 arcmin per panel\n'
                 'Cyan circle marks the 45-arcsec LoDeSS restoring beam at the CSS position')
    fig.tight_layout()
    fig.savefig(out/'j115948.png',dpi=150)
    plt.close(fig)
    pd.DataFrame(results).to_csv(out/'j115948_apertures.csv',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('out',type=Path)
    args=p.parse_args()
    inspect(args.root,args.out)
