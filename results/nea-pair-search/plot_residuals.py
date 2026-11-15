"""Plot observed radial residuals from the fully withheld dates and stations."""
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

plt.rcParams['svg.hashsalt'] = 'Open-Space-VE138'
plt.rcParams['svg.fonttype'] = 'none'
root = Path(__file__).resolve().parent
values = json.loads((root / 'pressure/fits/VE138.json').read_text())
with (root / 'orbits/fits/VE138.json').open() as stream:
    gravity = json.load(stream)
def rms(record):
    return record['radial_rms_arcsec']
categories = ['All 110 positions', 'Two stations held out (14)', 'Last 2007 night held out (3)',
              'Last 2026 night held out (2)', 'Later 2026 positions held out (59)']
grav = [rms(gravity['full']['training']), rms(gravity['hold']['withheld'])]
accel = [rms(values['all/A12']['training']), rms(values['sites/A12']['withheld'])]
for key in ('omit_last_old', 'omit_last_new', 'first_new8'):
    grav.append(rms(values[f'{key}/grav']['withheld']))
    accel.append(rms(values[f'{key}/A12']['withheld']))

fig, ax = plt.subplots(figsize=(8.2,3.3),layout='constrained')
y = np.arange(len(categories))
for yy, g, a in zip(y, grav, accel):
    ax.plot([g,a],[yy,yy], color='#adb8be', lw=2, zorder=1)
ax.scatter(grav,y,color='#3b5971',s=48,label='Gravity',zorder=2)
ax.scatter(accel,y,color='#bf5b27',marker='D',s=39,label='Gravity + acceleration',zorder=3)
ax.set_yticks(y, categories); ax.invert_yaxis()
ax.set_xlabel('Radial RMS positional error (arcseconds)')
ax.set_xlim(0,6.4); ax.xaxis.set_major_locator(MultipleLocator(1))
ax.grid(axis='x',color='#e0e5e8');ax.set_axisbelow(True)
ax.spines[['top','right','left']].set_visible(False)
ax.tick_params(axis='y',length=0,pad=9)
ax.legend(loc='upper right', frameon=False)
fig.savefig(root/'residuals.svg',metadata={'Date':None})
