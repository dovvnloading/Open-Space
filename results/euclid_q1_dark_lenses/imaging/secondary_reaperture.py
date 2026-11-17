"""Fit band-specific amplitudes at the independent mass geometry and recompute its aperture M/L."""
import argparse
import csv
import sys
import numpy as np
from scipy.ndimage import maximum_filter
from scipy.optimize import nnls
from scipy.signal import fftconvolve
from cutouts import CutoutArchive,image_extensions
from measure import HERE,blank_measurement,ab,table
from secondary_fit import ID,profiles
sys.path.insert(0,str(HERE.parent/'inference'))
from infer import mass_floor,luminosity

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--index-cache");args=p.parse_args()
    rec=next(r for r in csv.DictReader((HERE/'secondary_models.csv').open()) if r['model']=='SIS+shear')
    t=np.fromstring(rec['params'],sep=" ")
    a=CutoutArchive("lens",args.index_cache);d=a.members(ID,[ID+".fits","result/mge_lens_light.fits"])
    raw=image_extensions(d[ID+".fits"]);published_light=image_extensions(d["result/mge_lens_light.fits"])
    sun=np.genfromtxt(HERE.parent/"inference"/"sun_shifted_band_ab.csv",delimiter=",",names=True)
    yy,xx=np.indices((50,50));xx=xx-24.5;yy=yy-24.5
    crop=(slice(125,175),slice(125,175));ph,knots=[],[]
    for band in ['VIS','NIR_Y','NIR_J','NIR_H']:
        h,im=raw[band+'_FLUX'];err=raw[band+'_RMS'][1];q=raw[band+'_PSF'][1].astype(float);q/=q.sum()
        fg,src=0.,0.
        for dx in (-.25,.25):
            for dy in (-.25,.25):
                g,s=profiles(t,'SIS+shear',xx+dx,yy+dy);fg+=g/4;src+=s/4
        f=fftconvolve(fg,q,mode='same');s=fftconvolve(src,q,mode='same')
        good=np.isfinite(im)&np.isfinite(err)&(err>0)&(err<1e15)
        mask=(np.hypot(yy,xx)<21)&good[crop];design=np.column_stack((f[mask],s[mask]));w=1/err[crop][mask]
        data=im[crop][mask];centre=np.average(design,axis=0,weights=w*w);dmean=np.average(data,weights=w*w)
        amps,_=nnls((design-centre)*w[:,None],(data-dmean)*w);f*=amps[0];s*=amps[1]
        y,x=149.5+t[7],149.5+t[6]
        def measure(py,px,radius):
            dy,dx,bg,sd,_=blank_measurement(im,err,good,py,px,radius)
            ind=(round(py)+dy,round(px)+dx);loc=(ind[0]-125,ind[1]-125)
            obs=im[ind].sum()-bg;noise=max(sd,np.sqrt(np.sum(err[ind].astype(float)**2)))
            return obs,noise,f[loc].sum(),s[loc].sum()
        obs,sigma,fgf,sf=measure(y,x,t[8]);net=obs-sf;zp=float(h['MAGZERO']);zl=.17
        floor=mass_floor(t[8]*.1,zl)
        b='vis' if band=='VIS' else band[-1].lower()
        ph.append(dict(band=band,theta_arcsec=t[8]*.1,zl_assumed=zl,observed_counts=obs,source_counts=sf,
                       subtracted_counts=net,subtracted_snr=net/sigma,mag_E=ab(zp,net),mass_beta1=floor,
                       ml_beta1=floor/luminosity(ab(zp,net),zl,sun,b),ml_beta1_all_light=floor/luminosity(ab(zp,obs),zl,sun,b)))
        # Positions come solely from the VIS source model. The comparison in
        # the same data is a fit diagnostic, not a held-out counterimage test.
        if band=='VIS':
            rad=np.hypot(xx-t[6],yy-t[7])
            peaks=sorted(zip(*np.where((s==maximum_filter(s,5))&(rad>2)&(rad<20))),key=lambda p:-s[p])
            peaks=[pt for pt in peaks if s[pt]>.1*s.max()]
        for j,(iy,ix) in enumerate(peaks,1):
            obs,noise,mf,ms=measure(125+iy,125+ix,2.3)
            dy,dx,_,_,_=blank_measurement(im,err,good,125+iy,125+ix,2.3)
            pub=published_light[band][1][iy+50+dy,ix+50+dx].sum()
            knots.append(dict(band=band,peak=j,published_MGE_residual_snr=(obs-pub)/noise,dx_from_mass_arcsec=(xx[iy,ix]-t[6])*.1,dy_from_mass_arcsec=(yy[iy,ix]-t[7])*.1,
                              source_model_snr=ms/noise,foreground_subtracted_snr=(obs-mf)/noise))
    table(HERE/'secondary_reaperture.csv',ph);table(HERE/'secondary_refit_peaks.csv',knots)

if __name__=='__main__':main()
