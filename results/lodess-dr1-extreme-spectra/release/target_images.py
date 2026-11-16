"""Retrieve LoLSS/LoTSS images and compare 23, 54, 144 MHz for J112022+502936.

Usage: python target_images.py DATADIR
Uses the astrometric position of LOL1J112022.1+502936 in LoLSS DR1.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter

from acquire_lodess import retrieve, digest
from photometry import RadioImage
from resolution_check import aperture
from local_controls import nearby

def analyze(root):
    out=root/'j112022'
    out.mkdir(exist_ok=True)
    sky=SkyCoord(170.092166*u.deg,50.493238*u.deg)
    pos=sky.to_string(style='hmsdms',sep=':',precision=4)
    urls={'LoLSS-15':'https://lofar-surveys.org/public/lolss/mosaic-i/p170+50-mosaic.fits'}
    for key,path in [('LoTSS-6','dr3-cutout'),('LoTSS-20','dr3-low-cutout')]:
        urls[key]='https://lofar-surveys.org/'+path+'.fits?'+urlencode({'pos':pos,'size':12})
    def get(pair):
        k,url=pair
        f=out/(k+'.fits')
        retrieve(url,f)
        return dict(file=f.name,bytes=f.stat().st_size,sha256=digest(f),url=url)
    with ThreadPoolExecutor(max_workers=3) as pool:
        records=list(pool.map(get,urls.items()))
    pd.DataFrame(records).to_csv(out/'manifest.csv',index=False)
    files=[('LoDeSS',root/'images/P171+49-I.fits'),('LoLSS',out/'LoLSS-15.fits'),('LoTSS',out/'LoTSS-20.fits'),('LoTSS',out/'LoTSS-6.fits')]
    fig,axes=plt.subplots(1,4,figsize=(12.5,3.7)); res=[]
    for (survey,f),ax in zip(files,axes):
        with fits.open(f,memmap=True) as h:
            hdr=h[0].header
            beam=3600*hdr['BMAJ']; scale=3600*abs(hdr['CDELT1'])
            w=WCS(hdr).celestial
            s=Cutout2D(np.squeeze(h[0].data),sky,(7*u.arcmin,7*u.arcmin),wcs=w)
            arr=s.data.astype(float)
            x,y=s.wcs.world_to_pixel(sky)
            flux,bg,noise=aperture(arr,x,y,scale,beam)
            common=arr
            if beam<44.99:
                sig=np.sqrt(45.**2-beam**2)/(np.sqrt(8*np.log(2))*scale)
                common=gaussian_filter(arr-bg,sig,mode='nearest')*(45./beam)**2
            matched,_,_=aperture(common,x,y,scale,45)
            res.append(dict(survey=survey,beam_arcsec=beam,aperture_radius_arcsec=60,flux_jy=flux,matched_jy=matched,local_rms_jy_beam=noise))
            hdr_out=s.wcs.to_header()
            for k in ['BMAJ','BMIN','BPA','BUNIT','DATE-OBS','RESTFRQ','RESTFREQ','CRVAL3']:
                if k in hdr and k!='CRVAL3':hdr_out[k]=hdr[k]
            fits.writeto(out/(survey+f'-{beam:.0f}-stamp.fits'),s.data,hdr_out,overwrite=True)
            view=Cutout2D(arr-bg,(x,y),(100*u.arcsec,100*u.arcsec),wcs=s.wcs)
            ax.imshow(view.data,origin='lower',cmap='inferno',interpolation='nearest',vmin=-2*noise,vmax=max(5*noise,np.nanmax(view.data)*.8))
            ax.add_patch(Circle((.2*view.data.shape[1],.2*view.data.shape[0]),beam/(2*scale),fill=False,color='w',lw=1))
            ax.set_title(f'{survey} {int({"LoDeSS":23,"LoLSS":54,"LoTSS":144}[survey])} MHz\n{beam:.0f}″ beam',fontsize=10)
            ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle('J112022.1+502936 · 100″ field of view · individual brightness scales')
    fig.tight_layout();fig.savefig(out/'images.png',dpi=160);plt.close(fig)
    pd.DataFrame(res).to_csv(out/'apertures.csv',index=False)
    obs=[]
    for field in ['P171+49','P181+55','P151+49']:
        with RadioImage(root/'images'/(field+'-I.fits')) as im:
            h=im.hdus[0].header
            cen=SkyCoord(h['CRVAL1']*u.deg,h['CRVAL2']*u.deg)
            obs.append(dict(field=field,date_utc=h['DATE-OBS'],offaxis_deg=sky.separation(cen).deg,**im.measure(sky.ra.deg,sky.dec.deg)))
    pd.DataFrame(obs).to_csv(out/'decametre.csv',index=False)
    ctrl=nearby(pd.read_csv(root/'reference_sources.csv'),'P171+49',sky.ra.deg,sky.dec.deg)
    ctrl.to_csv(out/'controls.csv',index=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);analyze(p.parse_args().root)
