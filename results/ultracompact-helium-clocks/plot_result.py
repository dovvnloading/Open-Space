from pathlib import Path
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle
v=json.load(open('measurement.json'));b=pd.read_csv('blocks.csv');cross=pd.DataFrame(json.load(open('crossband.json'))['exposures_used'])
fig,axes=plt.subplots(1,3,figsize=(11.5,3.7),layout='constrained')
ax=axes[0]
for key,clr in [('pre_erupt','#adb8bd'),('early','#1d6880'),('late','#bc592f')]:
 d=b[b.window==key];ax.plot(d.btjd,np.hypot(d.sin,d.cos),'o',ms=3,c=clr)
ax.set(xlabel='TESS time (BJD − 2,457,000)',ylabel='6-hour fitted semi-amplitude (electrons/s)',title='The clock persists after the gap')
ax=axes[1];r=v['windows']['early'];im=ax.imshow(r['signed_image_e_s'],origin='lower',cmap='Blues')
for row in v['nearby_gaia']:
 if row['source_id'] in ['3458275544681382912','3458275643465151744','3458275544681383552']:
  star=row['source_id'][-4:];lab={'2912':'MASTER','1744':'Binary','3552':'Neighbor'}[star]
  ax.scatter(row['x'],row['y'],marker={'2912':'o','1744':'s','3552':'^'}[star],facecolors='none',edgecolors='#e05b2c',s=44,label=lab)
ax.set(xlim=(2.6,4.25),ylim=(2.35,4.5),xlabel='Detector column (pixels)',ylabel='Detector row (pixels)',title='Spatially resolved modulation')
ax.legend(fontsize=7,loc='upper right')
ax=axes[2]
for k,g in cross.groupby('night'):
 xx=g.tess_e_s-g.tess_e_s.mean();yy=g.magpsf.mean()-g.magpsf
 ax.errorbar(xx,yy,yerr=g.sigmapsf,fmt='o',ms=4,label=str(k))
x=np.array([-28,25]);ax.plot(x,x*json.load(open('crossband.json'))['slope_brightness_mag_per_e_s'],color='0.4',ls='--',lw=1)
ax.legend(title='UTC MJD',fontsize=7,title_fontsize=7);ax.set(xlabel='TESS modulation (electrons/s)',ylabel='ZTF brightening about nightly mean (mag)',title='Two observatories agree in phase')
for a in axes:a.tick_params(labelsize=8);a.title.set_fontsize(10);a.xaxis.label.set_fontsize(8);a.yaxis.label.set_fontsize(8)
fig.savefig('master-evidence.png',dpi=160)
