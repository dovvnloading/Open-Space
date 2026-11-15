"""Rebuild the folded photometry figure from the reported CSV measurements."""
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
with (root/'fold.csv').open() as handle:
    fold = list(csv.DictReader(handle))
with (root/'events.csv').open() as handle:
    events = list(csv.DictReader(handle))
fig, axes = plt.subplots(2,2,figsize=(9,5.8),layout='constrained')
for row, (name, ident, period) in enumerate([
    ('(29489) 1997 WQ','29489',12.353), ('(3177) Chillicothe','3177',18.560)]):
    records = [r for r in fold if r['mpc'] == ident]
    x = np.array([float(r['phase']) for r in records])
    y = np.array([float(r['median_mag']) for r in records])
    y -= np.median(y[abs(x)>.27])
    ax = axes[row,0]
    ax.axhline(0,color='#888888',lw=.7)
    ax.plot(x,y,'.-',color='#084B6B',lw=1.3,ms=5)
    ax.invert_yaxis()
    ax.set(ylabel='Relative magnitude',xlabel='Phase from dip center',title=f'{name}  ·  {period:.3f} h')
    visits = [r for r in events if r['mpc'] == ident]
    ax = axes[row,1]
    ax.axhline(0,color='#888888',lw=.7)
    for parity, color in [(0,'#084B6B'),(1,'#BE542B')]:
        xs=np.array([int(r['cycle']) for r in visits if int(r['cycle'])%2==parity])
        ys=np.array([float(r['delta_mag']) for r in visits if int(r['cycle'])%2==parity])
        ax.plot(xs,ys,'o',ms=3.5,color=color,label=['Even','Odd'][parity])
    ax.set(xlabel='Successive dip number',ylabel='Depth (mag)',title='Depth of individual dimmings')
    if row==0:
        ax.legend(frameon=False)
for ax in axes.flat:
    ax.spines[['top','right']].set_visible(False)
fig.savefig(root/'dimmings.png',dpi=160,facecolor='white')
