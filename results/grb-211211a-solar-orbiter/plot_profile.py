"""Plot the two independently recorded profiles at the predicted light-time offset."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent
data = pd.read_csv(root/'profiles.csv',dtype={'gbm_id':str})
p = data[data.gbm_id=='211211549']
t = p.seconds_from_predicted_arrival.to_numpy()
fig, (ax, bx) = plt.subplots(2, 1, figsize=(8,4.4), sharex=True, layout='constrained')
fig.suptitle('The hard pulse in STIX flare entry 2112111310',fontsize=15,x=.08,ha='left')
series = [(ax,p.stix_50_84,'#9a4300','Solar Orbiter / STIX','50–84 keV counts per 4 s'),
          (bx,p.gbm_channels_2_4_rate,'#1f6284','Fermi / GBM','Channels 2–4, summed rate')]
for axis, y, color, name, ylabel in series:
    axis.stairs(y,np.append(t-2,t[-1]+2),baseline=None,color=color,linewidth=2)
    axis.set_ylabel(ylabel,fontsize=9)
    axis.text(.98,.84,name,ha='right',transform=axis.transAxes,color=color,weight='bold')
    axis.spines[['top','right']].set_visible(False)
    axis.axvline(0,color='#808080',lw=1,ls=':')
    axis.grid(axis='y',alpha=.15)
ax.axhline(447,color='#9a4300',alpha=.4,lw=1,ls='--')
bx.set_xlabel('Seconds after Fermi trigger, with STIX shifted by the predicted 30.70 s')
fig.savefig(root/'profiles.png',dpi=150)
