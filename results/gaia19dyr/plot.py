"""Render the target measurements and the cross-instrument phase check."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent
d = pd.read_csv(root / 'nightly.csv')
d = d[d.star == 'target'].copy()
d['date'] = pd.to_datetime(d.mjd_night+2400000.5, origin='julian', unit='D')
fit = json.loads((root/'measurements.json').read_text())['fit']
fig, (ax, bx) = plt.subplots(2, 1, figsize=(8.4, 6.2), layout='constrained')
for band, label, color, size in [('zr','ZTF r','#405d9c',5),('G','Gaia G','#a73d32',14)]:
    q = d[d['filter']==band]
    ax.scatter(q.date,q.delta_mag,s=size,color=color,alpha=.7,label=label)
    h = q[q.date.dt.year.between(2019,2022)]
    phase=((h.jd-2459000)/fit['period_d']-fit['phase_cycle']+.5)%1-.5
    bx.scatter(phase,h.delta_mag,s=size,color=color,alpha=.6,label=label)
ax.set(ylabel='Relative magnitude',title='Gaia19dyr across the onset and decline of deep dips')
bx.set(xlabel='Phase of the 8.951-day fit to ZTF r (cycles)',ylabel='Relative magnitude',
       title='ZTF r and independent Gaia visits during 2019–2022')
bx.axvspan(-.16,.16,color='#eed58f',alpha=.32,zorder=-2)
for axis in (ax,bx):
    axis.invert_yaxis()
    axis.legend(frameon=False,loc='lower left',ncol=2)
    axis.spines[['right','top']].set_visible(False)
fig.savefig(root/'gaia19dyr.png',dpi=150)
