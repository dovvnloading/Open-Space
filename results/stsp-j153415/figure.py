"""Render the measured changes in separate DESI cameras."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from analyze import C, DATA, HERE, LINES, VISITS

fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
with np.load(DATA/'spectra.npz') as spec:
    for ax, name, title in zip(axes, ['Halpha','Hbeta','NaI'], ['Hα emission', 'Hβ emission', 'Na I absorption']):
        band, lines, _ = LINES[name]
        for exposure, label, color in zip(VISITS, ['2022 March 21', '2022 April 26'], ['#1c5b84','#b3481a']):
            key = f'{exposure}_{band}_'
            w = spec[key+'wave']
            flux = spec[key+'flux']
            usable = (spec[key+'mask']==0) & (spec[key+'ivar']>0)
            x = w-lines[0]
            select = usable & (x>-15) & (x<(27 if name=='NaI' else 15))
            sides = usable & ((abs(x+16)<3) | (abs(x-(lines[-1]-lines[0])-16)<3))
            ax.plot(x[select], flux[select]/np.median(flux[sides]), '.-', ms=3,
                    lw=1.5, label=label, color=color)
        for lam in lines:
            ax.axvline(lam-lines[0], color='#737373', lw=.8, linestyle=':')
        ax.set_title(title, loc='left', fontweight='bold')
        ax.set_xlabel(f'Wavelength minus {lines[0]:.2f} Å')
        ax.spines[['right','top']].set_visible(False)
axes[0].set_ylabel('Flux / local continuum')
axes[0].legend(frameon=False)
fig.savefig(HERE/'line-motion.png', dpi=160)
