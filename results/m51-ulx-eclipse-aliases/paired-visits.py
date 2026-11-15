"""Regenerate the source-resolved XMM rate plot from the measured windows."""
from pathlib import Path
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd

ROOT=Path(__file__).resolve().parent
initial=pd.read_csv(ROOT/'xmm'/'xmm-transition-checks.csv')
initial=initial[(initial['variation']=='base')&(initial.detector=='PN')].iloc[0]
next_visit=pd.read_csv(ROOT/'xmm'/'xmm-2018-recurrence.csv')
next_visit=next_visit[next_visit.detector=='PN'].set_index('window')
fig,axes=plt.subplots(1,2,sharey=True,figsize=(7.8,3.6),constrained_layout=True)
labels=['Measured rise','Predicted return 12.595 days later']
pairs=[([initial.ulx1_pre_rate,initial.ulx1_post_rate],[initial.ulx1_pre_error,initial.ulx1_post_error]),
       (next_visit.loc[['before','after'],'ulx1_rate'],next_visit.loc[['before','after'],'ulx1_error'])]
for ax,(rates,err),title in zip(axes,pairs,labels):
    ax.errorbar([0,1],rates,yerr=err,fmt='o',capsize=4,color='#8d3224',markersize=7)
    ax.set(xticks=[0,1],xticklabels=['Before','After'],xlim=(-0.55,1.55),ylim=(0,20),title=title)
    ax.grid(axis='y',alpha=.2)
    ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('M51 ULX-1 · pn counts per ks')
fig.savefig(ROOT/'xmm-repeat.png',dpi=145)
