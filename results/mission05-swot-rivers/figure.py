"""Regenerate the figures from the checked-in measurement tables."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import theilslopes

P = Path(__file__).resolve().parent
primary = pd.read_csv(P/'congo/paired-reach.csv')
held = pd.read_csv(P/'congo/holdout/visits.csv')
sat = primary.dropna(subset=['chenal_wse'])
test = held[held.period.eq('holdout')].dropna(subset=['chenal'])
fig, (ax, bx) = plt.subplots(1,2,figsize=(9.6,4.4),layout='constrained')
ink,blue,red='#192d37','#16757d','#b44a28'
for axis in [ax,bx]:
    axis.spines[['top','right']].set_visible(False)
    axis.grid(axis='y',color='#e0e5e5',lw=.8)
    axis.set_axisbelow(True)
for _,r in primary[primary.date.isin(['2024-01-08','2024-07-13'])].iterrows():
    high = r.date=='2024-01-08'
    ax.plot([0,1],[r.wse_upper,r.wse_lower],marker='o',lw=2.3,
            color=blue if high else red,label='8 Jan 2024 · flood' if high else '13 Jul 2024 · low water')
ax.set(xticks=[0,1],xticklabels=['Upstream reach','Downstream reach'],ylabel='SWOT water-surface elevation (m)',
       title='Water level rose farther downstream',xlim=(-.15,1.15),ylim=(243,275))
ax.legend(frameon=False,loc='lower left',fontsize=9)
slope,intercept,*_=theilslopes(sat.drop_m,sat.chenal_wse)
xx=np.linspace(sat.chenal_wse.min(),sat.chenal_wse.max())
bx.plot(xx,intercept+slope*xx,lw=1.4,color=ink,label='Fit to the first 17 visits')
bx.scatter(sat.chenal_wse,sat.drop_m,color=blue,alpha=.86,label='RiverSP 2.0 · 17 visits',s=34)
bx.scatter(test.chenal,test.drop_m,marker='D',facecolors='none',edgecolors=red,lw=1.5,label='Later RiverSP D · 7 visits',s=37)
bx.set(xlabel='Water elevation at the upstream Chenal reach (m)',ylabel='Outlet elevation drop (m)',
       title='The relationship recurred on later visits')
bx.legend(frameon=False,fontsize=8,loc='upper right')
for a in (ax,bx):
    a.title.set_fontsize(12)
fig.savefig(P/'congo-outlet.png',dpi=175,facecolor='white')
plt.close(fig)

profile=pd.read_csv(P/'congo/profile/profile.csv')
bins=pd.read_csv(P/'congo/profile/bins.csv').query('n>=3').copy()
fig,ax=plt.subplots(figsize=(8.8,4.0),layout='constrained')
for run,color,label in [('target',blue,'Jan–Jul 2024 · pass 111'),('check',red,'Dec 2023–Aug 2024 · pass 570')]:
    pts=profile[profile.run.eq(run)]
    seq=pts.node_idx+60*pts.reach_id_high.eq(13217000091)
    ax.scatter(seq,pts.rise_m,alpha=.28,s=11,color=color)
    group=bins[bins.run.eq(run)]
    centers=group.bin_start+2+60*group.reach_id_high.eq(13217000091)
    ax.plot(centers,group.rise_m,color=color,lw=1.9,label=label)
ax.axvspan(51,65,color='#414c54',alpha=.08)
ax.axvline(60.5,color='#7b8386',ls='--',lw=1)
ax.set(ylabel='High-water minus low-water elevation (m)',
       xlabel='SWORD v16 node sequence, downstream to upstream (nominal spacing 200 m)',
       title='Most of the change in water-surface drop occurs near the reach boundary')
ax.spines[['top','right']].set_visible(False)
ax.set_axisbelow(True)
ax.grid(axis='y',color='#e0e5e5')
ax.legend(frameon=False,fontsize=9)
fig.savefig(P/'outlet-profile.png',dpi=175,facecolor='white')
plt.close(fig)
