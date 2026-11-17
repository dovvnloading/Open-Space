"""Plot measured injection fractions using the deposited per-bin counts."""
import csv
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

HERE=Path(__file__).resolve().parent
rows=list(csv.DictReader((HERE/'tables/recovery.csv').open()))
fig,ax=plt.subplots(figsize=(6.4,3.65),layout='constrained')
colors={6:'#1c7884',10:'#a84d33'}
for length,label in [(6,'6 arcsec'),(10,'10 arcsec')]:
    selected=[r for r in rows if float(r['length_arcsec'])==length]
    x=[float(r['mu_ab_arcsec2']) for r in selected]
    y=[int(r['new_profile_match'])/int(r['attempts']) for r in selected]
    ax.plot(x,y,marker='o' if length==6 else 's',ms=6,color=colors[length],lw=1.8,label=label)
ax.set(xlim=(22.8,27.2),ylim=(-.04,1),xticks=range(23,28),xlabel='Mean observed VIS brightness within the FWHM (AB mag arcsec⁻²)',
       ylabel='New profile matches / injected',title='Straight 0.20-arcsec-wide wakes on real Euclid backgrounds')
ax.yaxis.set_major_formatter(PercentFormatter(1))
ax.grid(axis='y',color='#dedede',lw=.7)
ax.spines[['top','right']].set_visible(False)
ax.legend(frameon=False,title='Injected length')
fig.savefig(HERE/'recovery.png',dpi=165)
