"""Fit the positive point source on individual ZTF difference-image cutouts."""
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from astropy.coordinates import SkyCoord

# Pan-STARRS DR2 mean-object coordinate and 2025 TNS announcement position.
NUCLEUS = SkyCoord(350.42242651260784,26.48596952298109,unit='deg')
SN2025 = SkyCoord('23h21m41.388s','+26d29m11.23s')

def measure(image, wcs):
    image = image.astype(float)
    x0,y0 = wcs.world_to_pixel(NUCLEUS)
    y,x = np.mgrid[:image.shape[0],:image.shape[1]]
    ann = (np.hypot(x-x0,y-y0)>6)&(np.hypot(x-x0,y-y0)<17)
    bg = np.nanmedian(image[ann]); sigma=1.4826*np.nanmedian(abs(image[ann]-bg))
    use = (np.hypot(x-x0,y-y0)<5.5)&np.isfinite(image)
    xx,yy=x[use]-x0,y[use]-y0
    target=image[use]-bg
    def residual(z):
        amp,dx,dy,sx,sy,offset,bx,by=z
        model=amp*np.exp(-.5*(((xx-dx)/sx)**2+((yy-dy)/sy)**2))+offset+bx*xx+by*yy
        return (target-model)/sigma
    z=least_squares(residual,[max(1,target.max()),0,0,1.2,1.2,0,0,0],
                    bounds=([0,-3,-3,.55,.55,-np.inf,-np.inf,-np.inf],
                            [np.inf,3,3,2.5,2.5,np.inf,np.inf,np.inf]),loss='soft_l1').x
    pos=wcs.pixel_to_world(x0+z[1],y0+z[2])
    da,dd=NUCLEUS.spherical_offsets_to(pos)
    return dict(dx_arcsec=da.arcsec,dy_arcsec=dd.arcsec,sep_nucleus=pos.separation(NUCLEUS).arcsec,
                sep_2025=pos.separation(SN2025).arcsec,peak_over_rms=z[0]/sigma,
                sigma_x_pix=z[3],sigma_y_pix=z[4])

if __name__=='__main__':
    root=Path(__file__).resolve().parent
    rows=[]
    for filename in sorted(root.glob('cutouts_*.npz')):
        with np.load(filename, allow_pickle=False) as packet:
            for key in sorted(k for k in packet.files if k.startswith('img')):
                candid=key[3:]
                header=fits.Header.fromstring(str(packet['wcs'+candid]),sep='\n')
                result=measure(packet[key],WCS(header))
                rows.append(dict(candid=candid,**result))
    out=pd.DataFrame(rows).sort_values('candid')
    out.to_csv(root/'centroids.csv',index=False,float_format='%.7f')
    print(out.median(numeric_only=True).round(3).to_dict())
