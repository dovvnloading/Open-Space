"""Plot the measured F606W light curve and the resolved detector artifact."""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from measure import TARGET

rows = [r for r in csv.DictReader(open('visits.csv')) if r['band'] == 'F606W']
visits = np.array([int(r['visit']) for r in rows]); t = np.array([float(r['mjd_start']) for r in rows]); flux = np.array([float(r['target_e_s']) for r in rows]);
files = sorted(Path('stamps').glob('*020_drc.fits.npz'))
images = np.array([np.load(p)['image'] for p in files]); diffs = images - np.nanmedian(images[:41], axis=0)
fig = plt.figure(figsize=(9, 5.2), layout='constrained'); grid = fig.add_gridspec(2, 4, height_ratios=[1.15, 1]); ax = fig.add_subplot(grid[0, :]);
ax.scatter(t-t[51], flux, s=14, color='#526071'); ax.axhline(0,color='#aeb6c2',linewidth=.7); ax.set(xlabel='Days from observed F606W maximum', ylabel='Difference flux (electrons/s)',title='M87 nova 66 · HST/WFC3 optical photometry');
for i in [42,51]:
    ax.scatter(t[i]-t[51],flux[i],s=27,color='#a13320'); ax.annotate(f'Visit {i+1}',(t[i]-t[51],flux[i]),xytext=(-45,8),textcoords='offset points',fontsize=8)
for i, visit in enumerate([42,43,46,52]):
    a=fig.add_subplot(grid[1,i]); a.imshow(diffs[visit-1][32:61,34:63],origin='lower',cmap='gray',norm=PowerNorm(.55,vmin=-.045,vmax=.27)); a.set_xticks([]); a.set_yticks([]); a.set_title(f'Visit {visit}',fontsize=10)
    a.add_patch(plt.Circle((TARGET[0]-34,TARGET[1]-32),3.5,fill=False,color='#e87732',lw=1));
    if visit==43: a.annotate('artifact',(15,8),xytext=(.4,.05),textcoords='axes fraction',color='#f6d775',fontsize=8,arrowprops={'arrowstyle':'->','color':'#f6d775','lw':.7})
fig.savefig('lightcurve.png',dpi=140)
