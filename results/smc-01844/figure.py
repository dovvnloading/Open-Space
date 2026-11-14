"""Render the independently measured pulsation profiles and TESS frequency checks."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from replay import HERE, EPOCH, fit, series

m=json.loads((HERE/'measurements.json').read_text())
f=m['frequency_per_d']
fig, axes=plt.subplots(1,2,figsize=(11,5.3),gridspec_kw={'width_ratios':[1.38,1]},layout='constrained')
ax,bx=axes
selected={'OGLE-IV I':('OGLE-IV I',0,'#2058a2'),'EROS sm0066k12361 R':('EROS-2 red',.038,'#aa562f'),
          'Mt John IC':('Mt John I',.076,'#187b6e')}
for label,t,y,e,n in series():
    if label not in selected:
        continue
    name,shift,color=selected[label]
    b,cost,loss,r=fit(t,y,e,f,n)
    phase=(t-EPOCH)*f%1
    mod=np.column_stack([np.cos(2*np.pi*phase),np.sin(2*np.pi*phase)])@b[-2:]
    resid=r+mod
    bins=np.floor(phase*14).astype(int)
    centers,med=[] ,[]
    for i in range(14):
        values=resid[bins==i]
        if len(values)>=4:
            centers.append((i+.5)/14)
            med.append(np.median(values))
    ph=np.linspace(0,1,250)
    line=np.column_stack([np.cos(2*np.pi*ph),np.sin(2*np.pi*ph)])@b[-2:]
    ax.plot(ph,shift+line,color=color,lw=1.5)
    ax.scatter(centers,shift+np.array(med),color=color,s=20,label=name,zorder=3)
ax.set(xlabel='Phase at 4.88575 hours',ylabel='Residual magnitude (offset by series)',
       title='Ground-based photometry',xlim=(-.02,1.02))
ax.invert_yaxis();ax.legend(frameon=False,fontsize=9,loc='lower right')
v=json.loads((HERE/'tess.json').read_text())
for i,(s,d) in enumerate(v.items()):
    bx.plot([d['top_frequency'],f],[i,i],color='#9cabbc',lw=1.5)
    bx.scatter(d['top_frequency'],i,color='#2058a2',zorder=3,s=36)
bx.axvline(f,color='#a65233',linestyle='--',label='OGLE-IV frequency')
bx.set(yticks=np.arange(len(v)),yticklabels=list(v),ylabel='TESS sector',xlabel='Frequency (cycles per day)',
       xlim=(4.895,4.925),title='Highest TESS peak in each 0.1–15 d⁻¹ scan')
bx.legend(frameon=False,fontsize=9);bx.grid(axis='y',alpha=.18);bx.invert_yaxis()
for a in axes:
    a.spines[['top','right']].set_visible(False)
fig.savefig(HERE/'oscillation.png',dpi=165)
