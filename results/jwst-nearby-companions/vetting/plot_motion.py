"""Plot Pan-STARRS offsets from the 2025 F444W position and the host motion."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time
from matplotlib.ticker import MaxNLocator
from motion import SKY, PM, REF_MJD, parallax_shift
HERE=Path(__file__).resolve().parent
data=list(csv.DictReader((HERE/'ps1_motion.csv').open()))
t=np.array([float(r['mjd']) for r in data])
xy=np.array([[(float(r['ra'])-SKY[0])*np.cos(np.deg2rad(SKY[1]))*3600,(float(r['dec'])-SKY[1])*3600] for r in data])
grid=np.linspace(t.min()-25,REF_MJD,220)
pred=np.outer((grid-REF_MJD)/365.25,PM)+np.array([parallax_shift(v)-parallax_shift(REF_MJD) for v in grid])
fig,axs=plt.subplots(2,1,sharex=True,figsize=(8,5),layout='constrained')
for k,ax in enumerate(axs):
    ax.axhline(0,c='#555555',linewidth=1,label='Stationary background')
    ax.plot(Time(grid,format='mjd').jyear,pred[:,k],c='#b24822',label='Source following L 788−34')
    ax.scatter(Time(t,format='mjd').jyear,xy[:,k],s=20,c='#18627a',zorder=4,label='Pan-STARRS detections')
    ax.scatter(Time(REF_MJD,format='mjd').jyear,0,marker='*',c='#222222',s=95,zorder=5,label='JWST F444W')
    ax.spines[['right','top']].set_visible(False);ax.yaxis.set_major_locator(MaxNLocator(4));ax.grid(axis='y',alpha=.17)
    ax.set_ylabel(['East offset (arcsec)','North offset (arcsec)'][k])
axs[0].set_title('L 788−34 field source relative to its October 2025 position',loc='left')
axs[1].set_xlabel('Observation year')
axs[0].legend(frameon=False,ncol=2,fontsize=8,loc='lower right')
fig.savefig(HERE/'motion.png',dpi=160)
