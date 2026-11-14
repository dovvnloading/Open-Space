"""Plot predictions based entirely on pre-2022 data and the independent TESS observations."""
import json
from datetime import datetime, timedelta
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from joint import HERE, HOURS, timed_prediction
results=json.loads((HERE/'fits.json').read_text())
tess=json.loads((HERE/'tess_measurements.json').read_text())
def decimal_year(jd):
    dt=datetime(2000,1,1,12)+timedelta(days=float(jd)-2451545)
    return dt.year+(dt-datetime(dt.year,1,1))/(datetime(dt.year+1,1,1)-datetime(dt.year,1,1))
date=np.linspace(2459580.5,2461041,180)
x=[decimal_year(h) for h in date]
colors={'kepler':'#006c5d','single_trough':'#a64a21','cubic':'#615e87'}
names={'kepler':'Precessing orbit','single_trough':'One nonrepeating trough','cubic':'Cubic trend'}
fig,ax=plt.subplots(figsize=(7.3,3.9),layout='constrained')
for k,cc in colors.items():
    pred=HOURS*timed_prediction(k,np.array(results['before_2022'][k]['parameters']),date)
    ax.plot(x,pred,label=names[k],lw=2.0,color=cc)
for j,k in enumerate(['2023','2025']):
    d=tess[k];x0=decimal_year(np.mean(d['hjd_by_orbit']))
    p=HOURS*d['phase'];ci=np.quantile(d['bootstrap'],[.025,.975])*HOURS
    ax.errorbar(x0,p,yerr=np.array([[p-ci[0]],[ci[1]-p]]),fmt='o',color='#17242f',ms=6,
                capsize=4,label='Withheld TESS measurement (95% interval)' if j==0 else None)
ax.set(xlabel='Year',ylabel='Primary-to-secondary interval (hours)',xlim=(2021.95,2026),ylim=(14.68,16.29))
ax.set_xticks([2022,2023,2024,2025,2026]);ax.grid(axis='y',color='#d9dcde',linewidth=.6)
ax.spines[['right','top']].set_visible(False);ax.legend(loc='upper left',frameon=False,fontsize=8.2)
fig.savefig(HERE/'chronological-test.png',dpi=180)
