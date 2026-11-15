"""Draw the catalogue trajectories and the member-level separation estimates."""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from replay import HERE,INPUTS

rows = json.loads((HERE/'encounters.json').read_text())['candidates']
order = ['UPK_90','ASCC_99','NGC_6633','ASCC_114','Alessi_170']
rows.sort(key=lambda r:order.index(r['first']))
fig,(left,right)=plt.subplots(1,2,figsize=(11,4.3),gridspec_kw={'width_ratios':[1,1.25]})
blue='#194c75';rust='#ad4f20'
lead=rows[0]
for key,ls,col in [('Hunt & Reffert','-',blue),('UCC','--',rust)]:
    if key=='Hunt & Reffert':
        t,d,v=lead['t_myr'],lead['pericentre_pc'],lead['speed_kms']
    else:
        from replay import phase_space,encounter
        c=pd.read_csv(INPUTS/'ucc_centres.csv').set_index('name')
        states=[]
        for name in ['upk90','upk93']:
            s=c.loc[name]
            states.append(phase_space(s.RA_ICRS,s.DE_ICRS,1000/(s.Plx+.017),s.pmRA,s.pmDE,s.Rv))
        t,d,_,v=encounter(*states)
    grid=np.linspace(0,6,400)
    sep=lambda time:np.hypot(d,(time-t)*v*1.022712165)
    left.plot(grid,sep(grid),ls,color=col,label=key)
    if key=='Hunt & Reffert':
        epochs=np.arange(0,7,2)
        left.scatter(epochs,sep(epochs),s=22,c=col,zorder=3,label='2 Myr sampling')
left.axhline(20,linewidth=.8,color='#777777',linestyle=':',label='20 pc search threshold')
left.set(title='UPK 90 and UPK 93',xlabel='Time from present (Myr)',ylabel='Centre separation (pc)',ylim=(0,92),xlim=(0,6))
left.legend(frameon=False,fontsize=8)
for i,row in enumerate(rows):
    lo,med,hi=row['member_bootstrap'][1]['pericentre_pc']
    right.scatter(row['jacobi_sum_pc'],i,marker='|',color=rust,s=160,label='Sum of Jacobi radii' if i==0 else None)
    right.errorbar(med,i,xerr=[[med-lo],[hi-med]],color=blue,fmt='o',capsize=3,label='UCC median and 16–84%' if i==0 else None)
right.set_yticks(range(len(rows)),[r['first'].replace('_',' ')+' / '+r['second'].replace('_',' ') for r in rows])
right.invert_yaxis()
right.set(title='Five shortlisted passages',xlabel='Closest separation (pc)',xlim=(0,50))
right.legend(frameon=False,fontsize=8,loc='lower right',bbox_to_anchor=(1,.15))
for ax in (left,right):
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='x',color='#eeeeee')
fig.tight_layout(w_pad=2.1)
fig.savefig(HERE/'cluster_encounters.png',dpi=155,bbox_inches='tight')
