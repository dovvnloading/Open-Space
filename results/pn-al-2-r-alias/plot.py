"""Draw the measured clock and its sampling replica. Run beside analyze.py."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from analyze import LONG, SHORT, ROOT, harmonic, sinusoid


t, y, _ = np.loadtxt(ROOT / 'inputs/ogle-I.dat').T
_, slow = sinusoid(t, y, [LONG])
phase = (t - t[0]) / SHORT % 1
grid = np.linspace(0, 1, 300)
H = np.column_stack([np.ones(len(grid)), np.cos(2*np.pi*grid), np.sin(2*np.pi*grid)])
fig, axes = plt.subplots(2, 1, figsize=(7.0, 6.9), layout='constrained')
for series, color, text, style in [
    (y, '#243B63', 'OGLE observations', '-'),
    (slow, '#C14C19', '501.418-day model sampled at OGLE times', '--'),
]:
    coefficients, _ = sinusoid(t, series, [SHORT])
    axes[0].plot(grid, H @ coefficients, linestyle=style, lw=2.5, color=color, label=text, zorder=3)
axes[0].scatter(phase, y, s=6, alpha=0.15, c='#526885', linewidths=0)
axes[0].set(xlabel='Phase at 0.9952767 days', ylabel='OGLE I magnitude',
            title='The sampling alias reproduces the apparent short cycle', xlim=(0, 1))
axes[0].invert_yaxis()
axes[0].legend(loc='lower right', fontsize=8, framealpha=0.94)
H_short = harmonic(t, SHORT)
groups = np.cumsum(np.r_[True, np.diff(t)>0.5])
x, obs = [], []
coeff, _ = sinusoid(t, y, [SHORT])
for g in np.unique(groups):
    keep = groups == g
    if keep.sum() > 1:
        pred = H_short[keep] @ coeff[1:]
        x.extend(pred - pred.mean())
        obs.extend(y[keep] - y[keep].mean())
axes[1].axline((0,0), slope=1, color='#C14C19', ls='--', label='If the standalone short fit were physical')
axes[1].axhline(0, color='0.6', lw=1)
axes[1].scatter(x, obs, s=12, c='#243B63', alpha=.65, linewidths=0)
axes[1].set(xlabel='Short-cycle prediction within a night (mag)', ylabel='Measured change within a night (mag)',
            title='Repeat exposures reject the predicted variations', xlim=(-.135,.135), ylim=(-.135,.135))
axes[1].legend(loc='upper left', fontsize=8)
fig.savefig(ROOT / 'alias.png', dpi=150)
