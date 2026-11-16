"""Plot the concurrent RiverSP D visits and July–January staff-gauge comparisons."""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
s = pd.read_csv(ROOT/'analysis/concurrent-stages.csv')
s = s[(s.version=='D') & (s.floor==.75)]
pairs = pd.read_csv(ROOT/'analysis/matched-limbs.csv')
staff = pd.read_csv(ROOT/'external/hydrometry/staff_stage.csv')
fig, axes = plt.subplots(1,2, figsize=(10.5,4.3),gridspec_kw={'width_ratios':[1.3,1]})
colors = ['#087b86','#ad5820']
for (season,g), color in zip(s[s.season.isin([2024,2025])].groupby('season'), colors):
    g=g.sort_values('time_mopti')
    label=f'{season}–{str(season+1)[2:]}'
    axes[0].plot(g.mopti_m,g.dire_m,'-o',ms=3,lw=1.25,alpha=.85,color=color,label=label)
    p=pairs[pairs.season==season].iloc[0]
    v=g[g.time_mopti.isin([p.rising,p.receding])]
    axes[0].scatter(v.mopti_m,v.dire_m,s=58,color=color,edgecolors='#182127',linewidths=1.0,zorder=3)
    axes[0].plot(v.mopti_m,v.dire_m,ls=':',color=color,lw=1)
    gauge_dates=('2024-07-21','2025-01-26') if season==2024 else ('2025-07-13','2026-01-19')
    z=staff[(staff.station=='Diré') & staff.start.isin(gauge_dates)].sort_values('start')
    assert len(z)==2
    x=season-2024
    axes[1].bar(x-.17,p.dire_change_m,color=color,width=.31,hatch=None)
    axes[1].bar(x+.17,z.stage_cm.diff().iloc[-1]/100,facecolor='white',edgecolor=color,linewidth=1.8,width=.31,hatch='///')
    axes[1].text(x-.17,p.dire_change_m+.12,f'{p.dire_change_m:.2f}',ha='center',fontsize=9)
    axes[1].text(x+.17,z.stage_cm.diff().iloc[-1]/100+.12,f'{z.stage_cm.diff().iloc[-1]/100:.2f}',ha='center',fontsize=9)
axes[0].set(xlabel='Water-surface elevation at Mopti (m)',ylabel='Water-surface elevation at Diré (m)',title='Concurrent SWOT visits on pass 363')
axes[0].legend(frameon=False,fontsize=9)
axes[1].set(ylabel='January minus July stage at Diré (m)',title='Independent staff-gauge check',xticks=[0,1],xticklabels=['2024–25','2025–26'],ylim=(0,5.6))
from matplotlib.patches import Patch
axes[1].legend(handles=[Patch(facecolor='#53646b',label='SWOT'),Patch(facecolor='white',edgecolor='#53646b',hatch='///',label='Staff gauge')],frameon=False,ncol=2,loc='upper right',fontsize=9)
for ax in axes:
    ax.spines[['right','top']].set_visible(False)
    ax.grid(axis='y',alpha=.15)
    ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(ROOT/'matched-stage.png',dpi=145,bbox_inches='tight')
