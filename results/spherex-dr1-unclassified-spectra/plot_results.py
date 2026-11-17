"""Render the strongest off-pole peak and paired recovery on the deposited measurements."""
from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT=Path(__file__).resolve().parent
source='J135745.42+295804.4'
s=pd.read_csv(ROOT/'wide/spectra.csv').query('source_id == @source and 4.49 < wave_um < 4.72').sort_values('wave_um')
qa=pd.read_csv(ROOT/'wide/coadd_checks.csv').query("source_id == @source and radius_arcsec == 19 and sky_inner == 36 and mask == 'strict'").iloc[0]
repeats=json.loads((ROOT/'visits'/'wide'/(source+'.json')).read_text())['groups']['qr3_all']['bands']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.axisbelow':True})
fig,axes=plt.subplots(1,2,figsize=(10.3,4.25),layout='constrained')
ax=axes[0]
co=s.flux/qa.local_continuum_mJy
ax.errorbar(s.wave_um,co,s.sigma/qa.local_continuum_mJy,fmt='o-',lw=1.35,ms=4,color='#335c92',label='DR1 cube')
w=np.array([repeats[b]['wave'] for b in ['blue','core','red']])
v=np.array([repeats[b]['median'] for b in ['blue','core','red']])
e=np.array([repeats[b]['err'] for b in ['blue','core','red']])
base=np.interp(w,w[[0,2]],v[[0,2]])
ax.errorbar(w,v/base,e/base,fmt='D',ms=6,capsize=3,color='#b45728',label='Later QR3 images')
ax.axhline(1,color='#555',linewidth=.7,ls=':')
ax.set(xlabel='Observed wavelength (µm)',ylabel='Flux / local continuum',title='The strongest DR1 peak',xlim=(4.49,4.713))
ax.set_ylim(0, max(2.9,np.nanmax(co)+.2));ax.legend(frameon=False)
ax=axes[1]
for name,folder,color,marker in [('Six-detector field','wide/full','#335c92','o'),('Two-band polar field','nep_two_band','#24806d','s')]:
    t=pd.read_csv(ROOT/folder/'recovery_by_shape.csv').groupby('amplitude_frac')[['n','recovered']].sum()
    ax.plot(t.index*100,t.recovered/t.n,label=name,color=color,marker=marker,lw=1.8)
ax.set(xlabel='Injected peak / continuum (%)',ylabel='Localized recovery',title='Features injected into extracted spectra',
       ylim=(0,1.04),xlim=(6,103),xticks=[10,25,50,100])
ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0));ax.legend(frameon=False,loc='lower right')
for ax in axes:ax.grid(axis='y',color='#d9dfe4',lw=.6)
fig.savefig(ROOT/'search_results.png',dpi=165,bbox_inches='tight',facecolor='white')
