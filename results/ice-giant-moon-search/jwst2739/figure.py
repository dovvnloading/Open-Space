"""Plot the Proteus control and the full-population injection measurements."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

here=Path(__file__).parent
with (here/'halo/proteus_profile.csv').open() as f:
    profile=list(csv.DictReader(f))
r=json.loads((here/'plants_broad/summary.json').read_text())['recovery']
fig,axes=plt.subplots(1,2,figsize=(9.2,3.45),layout='constrained')
a,b=axes
for band,label,color in [('NIRCAM/F140M/CLEAR','F140M','#156f92'),('NIRCAM/F210M/CLEAR','F210M','#ae572d')]:
    use=[q for q in profile if q['band']==band]
    a.plot([float(q['radius_px']) for q in use],[float(q['max_score']) for q in use],label=label,color=color,linewidth=1.8)
a.axhline(5,color='0.4',linestyle=':',linewidth=1)
a.axvline(24,color='0.3',linestyle='--',linewidth=1)
a.set(xlabel='Separation from Proteus (output pixels)',ylabel='Largest co-moving control score',title='Light measured around Proteus')
a.legend(frameon=False)
flux=np.array([q['flat_uJy'] for q in r]); mag=-2.5*np.log10(flux*1e-6/3631)
cond=np.array([q['conditional_fraction'] for q in r]); interval=np.array([q['conditional_binomial_95'] for q in r]).T
b.errorbar(mag,100*cond,yerr=100*np.maximum(0,np.vstack((cond-interval[0],interval[1]-cond))),fmt='o-',capsize=3,color='#156f92',label='110 supported orbits')
b.plot(mag,[100*q['all_fraction'] for q in r],'s--',color='#ae572d',label='All 200 orbits')
b.axhline(50,color='0.4',linestyle=':',linewidth=1)
b.set(xlabel='Injected AB magnitude in both filters',ylabel='Recovered orbit fraction (%)',title='Single-visit recovery on real images',ylim=(-4,104))
b.legend(frameon=False,loc='lower left')
for ax in axes:
    ax.spines[['top','right']].set_visible(False)
fig.savefig(here/'neptune_completeness.png',dpi=145)
