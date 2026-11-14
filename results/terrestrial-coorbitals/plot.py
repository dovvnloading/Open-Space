"""Plot selected nominal orbits and the massless-Mars control from replay.py."""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

data = Path(sys.argv[1])
fig, axes = plt.subplots(2, 1, figsize=(8.4, 6.8), layout='constrained')
for label, color, caption in [
    ('mars-forward', '#00656a', 'JPL fit'),
    ('mars-esa-forward', '#be562d', 'ESA fit'),
    ('massless-mars', '#777777', 'Mars gravity removed'),
]:
    run = np.load(data / f'{label}.npz')
    x = run['time'] - run['time'][0]
    y = np.degrees(np.unwrap(np.radians(run['angle'][:, 0])))
    axes[0].plot(x, y, color=color, lw=1.6, label=caption)
axes[0].set(title='2026 RF64 relative to Mars', xlabel='Years after the 2026 orbit fit',
            ylabel='Unwrapped angle (degrees)')
axes[0].axhline(-60, lw=.6, color='#555555', ls=':')
h = np.genfromtxt(Path(__file__).with_name('horizons.csv'), delimiter=',', names=True)
x = 2000 + (h['jd_tdb'] - 2451545) / 365.25
y = (h['xn6_lambda_deg'] - h['emb_lambda_deg'] + 180) % 360 - 180
sel = (x > 2014) & (x < 2199)
axes[1].plot(x[sel], y[sel], color='#39394d', lw=2.7, label='JPL Horizons')
for label, color, caption in [('earth-forward', '#00656a', 'JPL fit'),
                              ('earth-esa', '#be562d', 'ESA fit')]:
    run = np.load(data / f'{label}.npz')
    active = run['time'] < 2200
    axes[1].plot(run['time'][active], run['angle'][active, 0], color=color,
                 lw=1.1, ls='--', label=caption)
axes[1].set(title='2025 XN6 relative to Earth', xlabel='Year',
            ylabel='Relative mean longitude (degrees)')
axes[1].axhline(0, color='#777777', lw=.8)
for ax in axes:
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', alpha=.15)
    ax.legend(loc='best', frameon=False, fontsize=9)
fig.savefig(sys.argv[2], metadata={'Date': None})
