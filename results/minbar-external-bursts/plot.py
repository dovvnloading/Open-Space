"""Plot the measured profiles at the lags from analyze.py."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from analyze import grouped

root = Path(__file__).resolve().parent
results = json.loads((root/'measurements.json').read_text())['waveforms']
xte = grouped('rxte_profiles.csv')
wind = grouped('konus_profiles.csv')
sources = [(2958,'2S 0918−549','2002 Jan 27'),(3165,'XTE J1810−189','2004 Feb 28'),
           (3238,'EXO 0748−676','2004 Nov 25'),(8279,'4U 1323−62','2011 Dec 05'),
           (3234,'4U 1702−429','2004 Nov 03 · solar flare')]

fig, axes = plt.subplots(5,1,figsize=(8,8.3),sharex=True,layout='constrained')
for ax,(seq,target,date) in zip(axes,sources):
    stat=results[str(seq)]['all']
    lag=stat['best_lag_s']
    t=np.array([float(x['start_s']) for x in xte[seq]])
    c=np.array([int(x['all']) for x in xte[seq]])
    kwt=np.array([float(x['start_s']) for x in wind[seq]])
    kwc=np.array([int(x['counts']) for x in wind[seq]])
    edges=np.arange(46.)
    inside=(t-lag>=0)&(t-lag<45)
    pca=np.histogram(t[inside]-lag,edges,weights=c[inside])[0]
    kon=np.histogram(kwt,edges,weights=kwc)[0]/np.histogram(kwt,edges)[0]
    for y,color,label in [(kon,'#156774','Konus-Wind, 50–200 keV'),
                           (pca,'#aa4325','RXTE/PCA, all accepted channels')]:
        ax.step(edges[:-1]+.5,(y-y.mean())/y.std(),where='mid',color=color,
                linewidth=1.5,alpha=.88,label=label)
    ax.set_title(f'{target}   ·   MINBAR {seq}   ·   {date}',loc='left',fontsize=10)
    ax.text(.99,.88,f"r = {stat['best_r']:.3f}",transform=ax.transAxes,ha='right',fontsize=9)
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(labelsize=8)
axes[-1].set_xlabel('Seconds from the Konus-Wind trigger')
fig.supylabel('Count rate relative to each curve’s mean (standard deviations)',fontsize=10)
axes[0].legend(loc='upper right',bbox_to_anchor=(1.0,1.59),ncol=2,fontsize=8,frameon=False)
fig.savefig(root/'waveforms.png',dpi=150)
