"""Draw conditional apparent positions in equatorial coordinates, with equal angular scale."""
from pathlib import Path
from datetime import date,timedelta
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from propagate import family,epoch,sky

here=Path(__file__).resolve().parent
model=family(json.loads((here/'orbit_input.json').read_text()))
start=date(2023,1,1);end=date(2026,12,31);step=8
dates=[start+timedelta(days=i) for i in range(0,(end-start).days+1,step)]
fig,ax=plt.subplots(figsize=(8.5,6.3),layout='constrained');cols=['#acd6d9','#79b9ba','#4c979e','#226c76'];lines=[]
sample=np.random.default_rng(62019).choice(len(model[0]),28,replace=False)
samplemodel=tuple(x[sample] for x in model[:3])+(model[3],)
for d in dates:lines.append(np.array(sky(samplemodel,epoch(d.isoformat()))))
lines=np.array(lines)
for y in range(2023,2027):
 m=np.array([d.year==y for d in dates]);ax.plot(lines[m,0],lines[m,1],color=cols[y-2023],lw=.75,alpha=.46)
for ds,label,c in [('2024-12-27','27 Dec 2024 (geocentric)','#674a91'),('2026-11-17','17 Nov 2026','#9a4527')]:
 points=np.array(sky(model,epoch(ds))).T;h=ConvexHull(points);reg=points[h.vertices];ax.fill(reg[:,0],reg[:,1],color=c,alpha=.19);ax.plot(*np.vstack([reg,reg[0]]).T,color=c,lw=1.6,label=label)
ax.set_xlabel('Right ascension (degrees, J2000)');ax.set_ylabel('Declination (degrees, J2000)');ax.set_title('Where the IRAS–AKARI pair would appear, 2023–2026',loc='left',fontsize=15,pad=13);ax.invert_xaxis();ax.set_aspect(1/np.cos(np.radians(-49.6)));ax.grid(alpha=.16);ax.spines[['top','right']].set_visible(False);ax.legend(loc='upper right',frameon=False);ax.ticklabel_format(useOffset=False)
fig.savefig(here/'search_region.png',dpi=140,facecolor='white');plt.close(fig)
from PIL import Image
file=here/'search_region.png'
with Image.open(file) as picture: compact=picture.convert('RGB').quantize(colors=128)
compact.save(file,optimize=True)
