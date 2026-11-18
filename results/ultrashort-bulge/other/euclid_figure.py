"""Render the three Euclid fields from real, registered VIS science cutouts."""
import csv
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import map_coordinates

from euclid_scene import NAMES, OUT, coords, cutout


def main():
    rows=list(csv.DictReader((OUT/'euclid_scene_checks.csv').open()))
    sources=list(csv.DictReader((OUT/'euclid_scene_sources.csv').open()))
    fig,axes=plt.subplots(1,3,figsize=(9.4,3.6),layout='constrained',sharey=True)
    for ax,(event,sex) in zip(axes,NAMES.items()):
        ra,dec=coords(*sex)
        # A common north-up, east-left grid. Interpolated masked pixels stay NaN.
        x=np.linspace(2.4,-2.4,121);y=np.linspace(-2.4,2.4,121)
        east,north=np.meshgrid(x,y)
        planes=[]
        for r in rows:
            if r['event']!=event or r['status']!='ok':continue
            h,im,flag,w,left,down=cutout(r['url'],ra,dec,[])
            im[(flag & ~((1<<18)|(1<<24)))!=0]=np.nan
            px,py=w.all_world2pix(ra+east/(3600*np.cos(np.deg2rad(dec))),dec+north/3600,0)
            planes.append(map_coordinates(im,[py-down,px-left],order=1,mode='constant',cval=np.nan))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            median=np.nanmedian(planes,axis=0)
        low,high=np.nanpercentile(median,(5,99.7));stretch=(high-low)/25
        image=np.arcsinh(np.maximum(median-low,0)/stretch)
        ax.imshow(image,origin='lower',extent=(2.42,-2.42,-2.42,2.42),cmap=plt.get_cmap('gray').with_extremes(bad='#574d68'),vmin=0,vmax=np.arcsinh(25))
        ax.plot(0,0,ls='',marker='+',color='#ffcc4d',mew=1.2,ms=12)
        for s in sources:
            if s['event']==event and int(s['independent_dithers'])>=3:
                ax.plot(float(s['east_arcsec']),float(s['north_arcsec']),ls='',marker='o',mfc='none',mec='#48d3ff',mew=.8,ms=8)
        ax.set_title(event,fontsize=10)
        ax.set(xlim=(2.4,-2.4),ylim=(-2.4,2.4),xticks=(2,1,0,-1,-2),yticks=(-2,-1,0,1,2),xlabel='East offset (arcsec)')
        ax.tick_params(labelsize=8)
    axes[0].set_ylabel('North offset (arcsec)')
    fig.savefig(OUT/'euclid-vis-scenes.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    main()
