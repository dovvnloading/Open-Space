"""Build the three Kepler postage-stamp panels from the checked-in pixel measurements."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ids = ["2079577161336953472","2079577157039120128","2079577161336952960"]
fig, axes = plt.subplots(2,3,figsize=(10.8,7.3))
fig.subplots_adjust(left=.09,right=.98,bottom=.13,top=.87,wspace=.14,hspace=.17)
palette = plt.get_cmap("magma").copy()
palette.set_bad("#c7cbd0")
for row,prefix,label in [(0,"","Q16 · KOI-7865"),(1,"q10_","Q10 · KOI-7866")]:
    maps=pd.read_csv(HERE/(prefix+"pixel_maps.csv"))
    field=pd.read_csv(HERE/(prefix+"gaia_field.csv"),dtype={"gaia_id":str}).set_index('gaia_id')
    filename = "q10_localizations.csv" if prefix else "eclipse_localizations.csv"
    measurements=pd.read_csv(HERE/filename).query('subset=="all"')
    panel_axes = axes[row]
    coords = [(int(r.x),int(r.y)) for r in maps.itertuples() if r.optimal_aperture]
    for ax,column,title in zip(panel_axes,['median_flux_e_s','primary_loss_e_s','secondary_loss_e_s'],
                              ['Median flux','Primary eclipse loss','Secondary eclipse loss']):
        shape = (int(maps.y.max())+1,int(maps.x.max())+1)
        data = np.full(shape,np.nan)
        data[maps.y,maps.x] = maps[column]
        ax.imshow(data,origin='lower',interpolation='nearest',cmap=palette,
                  norm=colors.PowerNorm(gamma=.5,vmin=0,vmax=np.nanmax(data)))
        for x,y in coords:
            for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)]:
                if (x+dx,y+dy) not in coords:
                    xc,yc=x+.5*dx,y+.5*dy
                    ax.plot([xc-.5*dy,xc+.5*dy],[yc-.5*dx,yc+.5*dx],color='white',alpha=.55,lw=.9)
        for sid,symbol,col in zip(ids,['+','o','x'],['#59d9f9','#58e8a0','#b5d2ff']):
            r=field.loc[sid]
            ax.plot(r.x_pix,r.y_pix,symbol,markeredgewidth=1.7,ms=10,color=col,markerfacecolor='none')
        ax.set(title=title if row==0 else None,xticks=range(shape[1]),yticks=range(shape[0]),xlim=(-.5,shape[1]-.5),ylim=(-.5,shape[0]-.5))
        ax.tick_params(length=0,labelsize=9)
        if row==1: ax.set_xlabel('Stamp column',fontsize=9)
        if ax==panel_axes[0]:ax.set_ylabel(label+'\nStamp row',fontsize=10)
        for side in ax.spines.values():side.set_visible(False)
    for ax,fit in zip(panel_axes[1:],measurements.itertuples()):
        ax.plot(fit.prf_x,fit.prf_y,'.',ms=4,color='white')
items=[Line2D([],[],ls='',marker=s,ms=8,color=c,markerfacecolor='none',label=l) for s,c,l in [
    ('+','#17829a','KOI-7865 target'),('o','#177448','Gaia eclipsing binary'),('x','#3e60a7','KOI-7866 target')]]
items += [Line2D([],[],ls='',marker='.',color='#333333',label='Fitted eclipse position'),
          Line2D([],[],ls='-',color='#999999',label='Kepler aperture')]
fig.legend(handles=items,loc='lower center',ncol=5,frameon=False,fontsize=9)
fig.suptitle('Two Kepler targets share one eclipsing source',fontsize=15,weight='bold')
fig.savefig(HERE/"kepler_difference_maps.png",dpi=175,facecolor='white')
