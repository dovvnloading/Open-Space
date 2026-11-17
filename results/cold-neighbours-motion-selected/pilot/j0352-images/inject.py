"""At eight off-source locations inject the same moving Gaussian into native W2 pixels.

This tests the local power of the dipole statistic at catalogue W2=16.923,
conditional on a circular FWHM 6.4-arcsec PSF and the published motion.
"""
import csv
import math
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import map_coordinates
from check import AX,CACHE,COS,DEC,FWHM,HERE,MA,MD,RA,X,Y,amp,manifest

def main():
    rows=manifest(); off=[(x,y) for x in (-18,0,18) for y in (-18,0,18) if (x,y)!=(0,0)]
    arrays=[];raw_arrays=[]
    for r in rows:
        with fits.open(CACHE/(r["file"]+".gz")) as f:
            h=f[0].header;data=f[0].data.astype(float)
        w=WCS(h).celestial;v,u=np.indices(data.shape)
        sx,sy=proj_plane_pixel_scales(w)*3600
        dt=(r["mjd"]-57170)/365.25
        amp_flux=10**(.4*(h["MAGZP"]-16.923))
        planted=data.copy()
        for ox,oy in off:
            cx,cy=w.world_to_pixel_values(RA+(ox+MA*dt)/(3600*COS),DEC+(oy+MD*dt)/3600)
            gaussian=np.exp(-4*np.log(2)*(((u-cx)*sx)**2+((v-cy)*sy)**2)/FWHM**2)
            planted+=gaussian*(amp_flux/gaussian.sum())
        xi,yi=w.world_to_pixel_values(RA+X/(3600*COS),DEC+Y/3600)
        for im,store in [(planted,arrays),(data,raw_arrays)]:
            a=map_coordinates(im,[yi,xi],order=1,cval=np.nan)
            store.append((a-np.nanmedian(a))*10**(.4*(20-h["MAGZP"])))
    out=[];summary=list(csv.DictReader((HERE/"test.csv").open()))
    for i,(early,late) in enumerate([([0,2],[18,20]),([1,3],[19,21])]):
        dts=(np.array([float(summary[i]["early_mjd"]),float(summary[i]["late_mjd"])])-57170)/365.25
        data=[]
        for images in [arrays,raw_arrays]:
            stacks=[np.nanmedian([a for a,r in zip(images,rows) if r["visit"] in visits],axis=0) for visits in [early,late]]
            data.append(stacks[1]-stacks[0])
        for ox,oy in off:
            a,b=[amp(diff,ox+MA*dts[0],oy+MD*dts[0],ox+MA*dts[1],oy+MD*dts[1]) for diff in data]
            out.append(dict(season=i,offset_east=ox,offset_north=oy,w2=16.923,injected_contrast=a,
                original_contrast=b,injected_increment=a-b,
                combined_contrast_over_mad=(a-float(summary[i]["off_median"]))/float(summary[i]["off_mad"])))
    with (HERE/"injections.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
    for i in range(2):
        vals=[r["combined_contrast_over_mad"] for r in out if r["season"]==i]
        print(i, np.median(vals), np.min(vals),np.max(vals),np.sum(np.array(vals)>=2),np.sum(np.array(vals)>=3))
if __name__=="__main__":main()
