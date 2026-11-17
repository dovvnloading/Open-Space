"""Show measured native visits and the DR1-trained ordinary PAH-band alternative."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p=Path(__file__).resolve().parent/'wide'
f=pd.read_csv(p/'pah_transfer.csv')
r=json.loads((p/'pah_transfer.json').read_text())
fig, axes=plt.subplots(1,2,figsize=(9.4,4.0),sharey=True,layout='constrained')
for ax, key, title in zip(axes,['dr1','qr3'],['DR1 · model fitted','QR3 · prediction from DR1']):
    g=f[f.release==key].sort_values('wave_um')
    unc=np.hypot(g.sigma,.03*abs(g.flux))
    ax.errorbar(g.wave_um,g.flux,yerr=unc,xerr=g.bandwidth_um/2,fmt='o',ms=5,color='#1f4570',
                lw=1.2,capsize=2,label='Native aperture')
    ax.scatter(g.wave_um,g.linear_plus_pah_prediction,s=58,marker='D',facecolor='none',
               edgecolor='#b84a21',lw=1.3,label='Continuum + PAH profile')
    grid=np.linspace(3.02,3.75,150)
    a,b=r['linear']['dr1_fit']
    ax.plot(grid, a+b*(grid-r['expected_um']),'--',color='#687482',lw=1.1,label='Continuum alone')
    ax.axvline(r['expected_um'],lw=.9,color='#b84a21',alpha=.55,ls=':')
    ax.set(xlim=(3.01,3.76),title=title,xlabel='Observed wavelength (µm)')
    ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('Background-subtracted aperture sum\n(MJy sr⁻¹ × pixel)')
axes[0].legend(loc='upper left',fontsize=8,frameon=False)
fig.suptitle('UGC 8856 NED01 · 3.29 µm at z = 0.03106',size=13)
fig.savefig(p/'pah_transfer.png',dpi=160)
