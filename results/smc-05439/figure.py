"""Fold the independent observations on the full alternating period."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from replay import HERE,EPOCH,design,templates,fit,training,records
import numpy as np

rows=list(records())
fs,b=training(rows)
fig,ax=plt.subplots(figsize=(7.6,5),layout='constrained')
display={'OGLE-4 I':('OGLE-IV I',0,'#22559a'),'OGLE-3 I':('OGLE-III I',.046,'#087a6d'),'EROS sm0044k13901 R':('EROS-2 red',.092,'#a95827')}
for name,t,y,e in rows:
    if name not in display:
        continue
    label,dy,col=display[name]
    a=np.column_stack((design(t,fs)[:,:19],templates(t,fs,b)[:,1:]))
    coeff,_=fit(a,y,e)
    resid=y-a[:,:19]@coeff[:19]
    ph=(t-EPOCH)*fs[1]%1
    g=np.linspace(0,1,301)
    model=templates(EPOCH+g/fs[1],fs,b)[:,1:]@coeff[-2:]
    ax.plot(g,dy+model,color=col)
    centers=[];med=[]
    for i in range(24):
        use=(ph>=i/24)&(ph<(i+1)/24)
        if sum(use)>=5:
            centers.append((i+.5)/24);med.append(np.median(resid[use]))
    ax.scatter(centers,dy+np.array(med),s=18,color=col,label=label)
ax.invert_yaxis();ax.legend(frameon=False);ax.set(xlabel='Phase at 1.165461 days',ylabel='Magnitude after subtracting the known binary (offset by series)',title='Two unequal dimmings toward OGLE-SMC-ECL-05439')
ax.spines[['top','right']].set_visible(False)
fig.savefig(HERE/'double-wave.png',dpi=165)
