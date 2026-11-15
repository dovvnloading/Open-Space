"""Plot the measured bands against the same-detector local continuum."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from analyze import BANDS, band_measurement

root = Path(__file__).resolve().parent
with (root/'spectrum.csv').open() as stream:
    table = list(csv.DictReader(stream))
wave = np.array([float(r['wavelength_um']) for r in table]); fig, panels=plt.subplots(1,3,figsize=(10.2,3.9))
targets=[('MN39','aperture_sum_mjy_sr','#005f73','o'),('J161036.07−512153.8','neighbor_aperture_sum_mjy_sr','#a34414','s')]
for axis,(label,band),title in zip(panels,BANDS.items(),['H₂O  ·  3.03 μm','CO₂  ·  4.26 μm','CO  ·  4.67 μm']):
    center,left,right=band
    for target,key,color,marker in targets:
        flux=np.array([float(row[key]) if row[key] else np.nan for row in table]); xp=[]; yp=[]
        for low,high in (left,right):
            use=(wave>=low-.0005)&(wave<=high+.0005)&np.isfinite(flux)
            xp.append(np.median(wave[use]));yp.append(np.log(np.median(flux[use])))
        valid=(wave>=left[0]-.0005)&(wave<=right[1]+.0005);continuum=np.exp(yp[0]+(wave-xp[0])*(yp[1]-yp[0])/(xp[1]-xp[0]));axis.plot(wave[valid],(flux/continuum)[valid],color=color,marker=marker,ms=4,lw=1.4,label=target)
    axis.axhline(1,color='#627177',lw=.8,ls=':');axis.set_xlabel('Wavelength (μm)');axis.set_title(title,loc='left',fontsize=12,fontweight='bold');axis.set_ylim(0,1.55);axis.spines[['top','right']].set_visible(False);axis.grid(axis='y',alpha=.18)
panels[0].set_ylabel('Signal / local continuum');panels[0].legend(frameon=False,fontsize=8,loc='upper right');fig.tight_layout();fig.savefig(root/'ice-bands.png',dpi=170,bbox_inches='tight');plt.close(fig)
