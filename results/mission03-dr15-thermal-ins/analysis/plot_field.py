"""Draw a CFHT r finding chart and the background-inclusive EPIC-pn source spectrum."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.patches import Circle
from scipy.ndimage import map_coordinates
from fit_pha import read, sigma, photon_shape

def figure(image,pha,fit_json,output):
    ra,dec=134.077254,38.119831
    data,h=fits.getdata(image,header=True)
    w=WCS(h).celestial
    axis=np.linspace(-13,13,251)
    west,north=np.meshgrid(axis,axis)
    x,y=w.world_to_pixel_values(ra-west/(3600*np.cos(np.deg2rad(dec))),dec+north/3600)
    stamp=map_coordinates(data,[y,x],order=1)
    fig,(sky,spec)=plt.subplots(1,2,figsize=(9.1,4.15),gridspec_kw={'width_ratios':[1,1.2]},layout='constrained')
    med=np.nanmedian(stamp);sig=1.4826*np.nanmedian(abs(stamp-med))
    sky.imshow(np.arcsinh((stamp-med)/sig),origin='lower',extent=(-13,13,-13,13),cmap='Greys',vmin=-1.5,vmax=3,interpolation='nearest')
    sky.add_patch(Circle((0,0),3.204,edgecolor='#287295',linewidth=1.8,fill=False,label='DR15 99% circle'))
    sky.add_patch(Circle((-.347,.949),3.204,edgecolor='#aa522b',linewidth=1.6,ls='--',fill=False,label='Gaia-shifted circle'))
    sky.set(xlabel='West offset (arcsec)',ylabel='North offset (arcsec)',title='CFHT MegaCam · r band')
    sky.legend(loc='lower left',fontsize=8,framealpha=.94)
    modelinfo=json.loads(fit_json.read_text())['models']
    en,mat,n,b,alpha,lo,hi,_,_=read(pha,.3)
    # Display summed source-region and background counts in broad energy bins.
    group_edges=np.array([.3,.4,.5,.6,.8,1.,1.5,2.,3.,5.])
    group=np.searchsorted(group_edges,(lo+hi)/2,side='right')-1
    def grouped(arr):return np.bincount(group,weights=arr,minlength=len(group_edges)-1)
    src=grouped(n);bg=grouped(alpha*b)
    center=(group_edges[1:]+group_edges[:-1])/2
    spec.errorbar(center,src,yerr=np.sqrt(src+.25),linestyle='none',marker='o',ms=4,color='#222222',label='Source aperture')
    spec.stairs(bg,group_edges,color='#868686',linewidth=1.4,label='Scaled background')
    for key,col,label in [('blackbody','#287295','Blackbody + background'),('powerlaw','#aa522b','Power law + background')]:
        info=modelinfo[key]
        photons=info['normalization']*photon_shape(en,key,info['shape'])*np.exp(-info['NH_1e21']*1e21*sigma(en))
        spec.stairs(grouped(mat@photons)+bg,group_edges,color=col,linewidth=1.6,label=label)
    spec.set(xscale='log',xlabel='Photon energy (keV)',ylabel='Counts per displayed bin',title='XMM-Newton · EPIC-pn')
    spec.set_xticks([.3,.5,1,2,3,5],['0.3','0.5','1','2','3','5']);spec.set_xlim(.29,5.1);spec.set_ylim(bottom=0)
    spec.legend(frameon=False,fontsize=8);spec.spines[['top','right']].set_visible(False)
    fig.savefig(output,dpi=165);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['image','pha','fit_json','output']:p.add_argument(key,type=Path)
    a=p.parse_args();figure(a.image,a.pha,a.fit_json,a.output)
