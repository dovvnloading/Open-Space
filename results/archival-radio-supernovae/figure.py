"""Plot the VLA pixels and fixed-position fits of the two strongest supernova candidates."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

root=Path(__file__).resolve().parent
tab=pd.read_csv(root/'beam-fits.csv').merge(pd.read_csv(root/'followup-visits.csv'),on=['name','epoch'])
fig=plt.figure(figsize=(11.8,5.6),layout='constrained')
grid=GridSpec(2,6,figure=fig,width_ratios=[1,1,1,1,.12,2.55])
for row,(name,cap) in enumerate([('SN2019njv',8.5),('SN2019bao',4.2)]):
    rows=tab[tab.name==name].sort_values('epoch')
    for col,r in enumerate(rows.itertuples(index=False)):
        h=fits.getheader(root/'inputs'/f'{name}-{r.epoch}.fits')
        im=fits.getdata(root/'inputs'/f'{name}-{r.epoch}.fits').squeeze()*1000
        w=WCS(h).celestial
        north,east=np.mgrid[-8:8.01:.2,8:-8.01:-.2]
        x,y=w.world_to_pixel_values(r.ra+east/3600/np.cos(np.deg2rad(r.dec)),r.dec+north/3600)
        img=map_coordinates(im,[y,x],order=1)
        ax=fig.add_subplot(grid[row,col])
        ax.imshow(img,origin='lower',extent=(-8,8,-8,8),cmap='magma',vmin=-.45,vmax=cap)
        # East points left and north points up. White ticks mark the optical coordinate.
        for sign in [-1,1]:
            ax.plot([sign*2.2,sign*3.6],[0,0],c='#72e0e4',lw=1.1)
            ax.plot([0,0],[sign*2.2,sign*3.6],c='#72e0e4',lw=1.1)
        ax.add_patch(Ellipse((-5.5,-5.5),r.beam_minor,r.beam_major,angle=r.beam_pa,fc='none',ec='white',lw=.7))
        ax.set_title(r.date_x[:10],fontsize=9)
        ax.set_xticks([]);ax.set_yticks([])
        if col==0:ax.set_ylabel(name.replace('SN','SN '),weight='bold',fontsize=12)
    ax=fig.add_subplot(grid[row,5])
    ax.axhline(0,c='.65',lw=.7);ax.axvline(0,c='.72',ls=':',lw=.9)
    ax.errorbar(rows.days/365.25,rows.forced,yerr=rows.scatter,fmt='o',color='#16566d',capsize=3,ms=5)
    ax.set_ylabel('3 GHz (mJy/beam)',fontsize=9)
    if row==1:ax.set_xlabel('Years from optical maximum')
    ax.spines[['top','right']].set_visible(False)
fig.suptitle('Late radio emission at two supernova positions',fontsize=16,weight='bold')
fig.savefig(root/'radio-history.png',dpi=160)
