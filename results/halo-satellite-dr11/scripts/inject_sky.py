"""Plot actual secondary injection positions for M_V=-2, r_h=30 pc, DM=19."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Ellipse, Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from off_control import secondary_grid
from search import FIELDS, CONTROLS, OUT, STEP


def run():
    trials = pd.read_csv(OUT/'off_control_injection_trials.csv.gz')
    selected = trials[(trials.mv_expected==-2)&(trials.rhalf_pc==30)&(trials.dm==19)]
    fig,axes = plt.subplots(1,2,figsize=(9.8,4.9))
    fig.subplots_adjust(left=.07,right=.99,bottom=.25,top=.81,wspace=.21)
    fig.suptitle('Where the injected faint systems were recovered',fontsize=15,y=.98)
    fig.text(.5,.89,'Expected $M_V=-2$  ·  projected half-light radius 30 pc  ·  distance 63 kpc',ha='center')
    for axis,(field,(ra0,ra1,dec0,dec1)) in zip(axes,FIELDS.items()):
        data = selected[selected.field==field]
        grid = secondary_grid(field)
        mask = grid.valid[0]
        xmax = grid.ra0+grid.nx*STEP/grid.cos
        ymax = grid.dec0+grid.ny*STEP
        axis.imshow(mask,origin='lower',interpolation='nearest',extent=[ra0,xmax,dec0,ymax],
                    cmap=ListedColormap(['white','#e4eaea']),vmin=0,vmax=1,aspect='auto')
        cr,cd,_ = CONTROLS[field]
        axis.add_patch(Ellipse((cr,cd),2*26/60/grid.cos,2*26/60,facecolor='none',edgecolor='#59606a',
                               lw=1,ls='--',zorder=2))
        axis.text(cr,cd,'Boötes I' if field=='bootes' else 'Coma Ber',ha='center',va='center',fontsize=8,color='#444e5b')
        for ok,marker,color in [(True,'o','#086955'),(False,'X','#aa4c24')]:
            group = data[data.recovered==ok]
            axis.scatter(group.ra,group.dec,marker=marker,s=37,c=color,edgecolors='white',linewidths=.6,zorder=3)
        xmid=(ra0+xmax)/2; ymid=(dec0+ymax)/2
        axis.axvline(xmid,lw=.65,color='#a5acb0');axis.axhline(ymid,lw=.65,color='#a5acb0')
        axis.set_xlim(ra1,ra0);axis.set_ylim(dec0,dec1)
        axis.set_aspect(1/grid.cos)
        axis.set_xlabel('Right ascension (deg)');axis.set_ylabel('Declination (deg)')
        name = 'Boötes' if field=='bootes' else 'Coma'
        axis.set_title(f'{name} field  ·  {data.recovered.sum()} / {len(data)} recovered',fontsize=11)
        axis.tick_params(labelsize=9)
    handles = [Line2D([],[],marker=m,ls='none',color=c,label=l) for m,c,l in
               [('o','#086955','Recovered'),('X','#aa4c24','Missed')]]
    handles.append(Patch(facecolor='#e4eaea',label='Eligible 3′ aperture centres'))
    handles.append(Line2D([],[],color='#59606a',ls='--',label='26′ known-system exclusion'))
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.08),ncol=4,frameon=False,fontsize=9)
    fig.text(.5,.02,'Six catalogue-level injections per quadrant. Eligible centroid area totals 7.366 deg² at 3′ and decreases for larger apertures.',ha='center',fontsize=8.7)
    with matplotlib.rc_context({'svg.fonttype':'none','svg.hashsalt':'dr11-secondary-sky-v1'}):
        fig.savefig(OUT/'off_control_sky_completeness.svg',metadata={'Date':None})
    plt.close(fig)


if __name__ == '__main__':
    run()
