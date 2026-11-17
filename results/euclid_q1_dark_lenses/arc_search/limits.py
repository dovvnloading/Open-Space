"""Redshift- and profile-dependent sensitivity for the morphology triggers.

All results assume the trial radius is an Einstein radius. Luminosity uppers refer
to observed (attenuated) light in a matched, blueshifted passband.
"""
import csv
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np
from scipy import ndimage as ndi,special
from search import HERE
from photometry import WIDTH,save_csv
sys.path.insert(0,str(HERE.parent/'inference'))
from infer import chi,mass_floor,luminosity

@lru_cache(None)
def fraction(band,re,n,r):
    # Oversampled PSF-convolved circular Sersic light in fixed core vs E. Regularize
    # the central pixel via an even half-pixel grid so the Sersic cusp is integrated.
    step=.02
    axis=(np.arange(-220,220)+.5)*step
    xx,yy=np.meshgrid(axis,axis)
    rr=np.hypot(xx,yy)
    bn=special.gammaincinv(2*n,.5)
    counts=np.exp(-bn*(rr/re)**(1/n))
    image=ndi.gaussian_filter(counts,.1*WIDTH[band]/(2.355*step))
    return image[rr<=.3].sum()/image[rr<=r].sum()

def main():
    data=list(csv.DictReader((HERE/'validation/photometry.csv').open()))
    solar=np.genfromtxt(HERE.parent/'inference/sun_shifted_band_ab.csv',names=True,delimiter=',')
    phot={(v['id'],v['band'],v['method']):v for v in data if v['point']=='center'}
    rows=[]
    for (name,b,ap),w in sorted(phot.items()):
        if ap!='native_E':continue
        r=float(w['radius_arcsec'])
        direct=float(w['flux_upper5_uJy'])
        core=float(phot[name,b,'native_core']['flux_upper5_uJy'])
        for n,re in [(0,0),(1,.2),(1,.8),(4,.2),(4,.8)]:
            # Report the two limits separately; no confidence statement for the
            # minimum of several profile-dependent bounds.
            f=fraction(b,re,n,r) if n else np.nan
            fupper=core/f if n else direct
            mag=23.900065622-2.5*np.log10(fupper)
            for zl in (.15,.3,.6,1.,1.4):
                zs=2.
                beta=1.-chi(zl)/chi(zs)
                mass=mass_floor(r,zl)/beta
                light=luminosity(mag,zl,solar,b.lower())
                rows.append(dict(id=name,band=b,theta_assumed_arcsec=r,zl_assumed=zl,zs_assumed=zs,
                    constraint='core_profile' if n else 'all_light_in_E',sersic_n=n,re_arcsec=re,
                    core_fraction_of_E=round(f,6),flux_upper5_uJy=round(fupper,6),m_E_lower_AB=round(mag,4),
                    projected_mass_Msun=round(mass),observed_L_E_upper_Lsun=round(light),
                    conditional_ML_floor=round(mass/light,5)))
    save_csv(HERE/'validation/conditional_limits.csv',rows)
if __name__=='__main__':main()
