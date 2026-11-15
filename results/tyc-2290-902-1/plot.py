"""Plot the spectra and independent photometric measurements used in the report."""
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
from analyze import C, ROOT, rows

result = json.loads((ROOT / 'results.json').read_text())
data = ROOT / 'data'
spectra = np.load(data / 'spectra.npz')
fig, axes = plt.subplots(1, 2, figsize=(11, 4.1), constrained_layout=True)
ax, phot = axes
w, s = spectra['684208204_COADD_B']
window = (w > 5160) & (w < 5300)
for row, color in zip(rows(data, 'controls.csv'), ['#095888', '#C06419', '#47826C']):
    wt, z = spectra[row['obsid'] + '_COADD_B']
    vel = np.arange(-320, 321, 2)
    values = [np.corrcoef(s[window], np.interp(w[window]/(1 + (v - float(row['rv_b1']))/C), wt, z))[0, 1] for v in vel]
    ax.plot(vel, values, color=color, lw=1.6, label=row['obsid'])
ax.set(title='Two stellar velocity peaks', xlabel='Trial velocity (km/s)', ylabel='Normalized cross-correlation')
ax.legend(title='Same-plate template', fontsize=8, frameon=False)
tess = np.genfromtxt(data / 'tess.csv', names=True, delimiter=',')
p = result['ztf_g_period']
for sector, color in zip([17, 57, 84], ['#095888', '#C06419', '#47826C']):
    a = tess[tess['sector'] == sector]
    x = ((a['btjd']-1700)/p) % 1
    vals = (a['flux']-np.median(a['flux']))*100
    centers = (np.arange(24)+.5)/24
    m = [np.median(vals[abs(x-c)<.5/24]) for c in centers]
    phot.plot(centers,m,'o-',color=color,lw=1.3,ms=3,label=f'TESS {sector}')
g = [r for r in rows(data,'ztf.csv') if r['oid']=='650115100004176']
y = np.array([float(r['mag']) for r in g]);t=np.array([float(r['hjd'])-2457000 for r in g]);use=np.array([r['catflags']=='0' and float(r['magerr'])<.06 for r in g])&(abs(y-np.median(y))<.15)
ls=LombScargle(t[use],y[use]);x=np.linspace(0,1,160);y=ls.model(1700+x*p,1/p);phot.plot(x,(10**(-.4*(y-y.mean()))-1)*100,'k--',lw=1.7,label='ZTF g fit')
phot.set(title='Three TESS visits folded on the ZTF period', xlabel='Phase at 0.535397 days', ylabel='Relative brightness (%)',xlim=(0,1));phot.legend(frameon=False,fontsize=8)
for a in axes:
    a.spines[['top','right']].set_visible(False)
fig.savefig(ROOT/'evidence.png',dpi=160)
