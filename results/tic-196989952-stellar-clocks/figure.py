from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from analyze import inputs, outside

root=Path(__file__).resolve().parent
result=json.loads((root/'measurements.json').read_text())
freq=result['tess_frequencies_cpd'][0]
tess,w,_=inputs()
fig,axes=plt.subplots(1,2,figsize=(10.2,3.8),layout='constrained')
colors=['#276a76','#b34229']
for j,(name,df) in enumerate([('TESS, sector 60',tess[tess.sector==60]),('SuperWASP, 2006–08',w)]):
    t=df.btjd.to_numpy();y=df.flux.to_numpy();y=y/np.median(y)
    use=outside(t)&(abs(y-1)<.09)
    t,y=t[use],y[use];phi=(t*freq)%1
    # Display independently folded shapes, aligning their fitted maxima for visualization.
    X=np.c_[np.ones(len(t)),np.cos(2*np.pi*phi),np.sin(2*np.pi*phi)]
    b=np.linalg.lstsq(X,y,rcond=None)[0]
    phase=(phi-np.arctan2(b[2],b[1])/(2*np.pi))%1
    bins=np.arange(20)/20+.025;fold=np.array([np.median(y[(phase>=k-.025)&(phase<k+.025)]) for k in bins])
    axes[0].plot(bins,(fold/np.mean(fold)-1)*100,'.-',ms=5,lw=1.4,label=name,color=colors[j])
axes[0].set(xlabel='Phase of the 0.84445-day signal',ylabel='Relative flux (%)',title='Out-of-eclipse modulation')
axes[0].axhline(0,lw=.6,color='.65');axes[0].legend(frameon=False,fontsize=9)
for j,row in enumerate(result['spacings']):
    x=(row['mean_btjd']+2457000-2451544.5)/365.25+2000
    y=row['separation_days']*24;lo,hi=np.array(row['bootstrap_p16_p84'])*24
    axes[1].errorbar(x,y,yerr=[[y-lo],[hi-y]],fmt='o',capsize=3,color=colors[1 if j==0 else 0],ms=6)
axes[1].set(xlabel='Mean observing year',ylabel='Primary-to-secondary separation (hours)',title='The eclipse separation grows by 3.24 hours')
axes[1].xaxis.set_major_locator(MaxNLocator(integer=True,nbins=5))
for a in axes:a.spines[['top','right']].set_visible(False)
fig.savefig(root/'clocks.png',dpi=160)
