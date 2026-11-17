"""Collect image-based projected phase offsets; compare observing epochs without fitting a drift law."""
import csv
from datetime import date,timedelta
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.dates as dates
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parent


def rows(path):
    with (ROOT/path).open() as f:return list(csv.DictReader(f))


def median(values):
    return float(np.median(list(values)))


def main():
    result=[]
    def add(mjd,source,phase,n,kind):
        utc=date(1858,11,17)+timedelta(days=int(mjd))
        result.append(dict(utc=utc.isoformat(),object=source,longitude_residual_deg=round(float(phase),4),n_positions=n,measurement=kind))
    for name in ('perdita','belinda'):
        data=rows(f'jwst2739_uranus/f140m_{name}_positions.csv')
        for visit in sorted(set(x['visit'] for x in data)):
            group=[r for r in data if r['visit']==visit]
            add(float(group[0]['mjd_utc']),name.title(),median(float(r['peak_longitude_deg']) for r in group),len(group),'F140M native trail fits')
    for key,group in json.loads((ROOT/'jwst2768/perdita/frame_longitudes.json').read_text()).items():
        add(60346,'Perdita' if key=='725' else 'Belinda',group['median_deg'],group['n'],'F182M and F210M across February 6–7')
    checks=json.loads((ROOT/'jwst6379_9482/psf_astrometry/ephemeris_checks.json').read_text())
    for mjd,k in ((60708,'first'),(60955,'late')):
        add(mjd,'Perdita',checks[k]['phase_diff_deg'],checks[k]['n'],'F150W2 native integrations')
    short=[float(r['phase_offset_deg']) for r in rows('nircam_supplement/f140m_catalog.csv') if r['id']=='725']
    long=[float(r['longitude_deg']) for r in rows('jwst8975/long_holdout/perdita_full_trails.csv') if r['accepted']=='True']
    add(61321,'Perdita',median(short+long),len(short)+len(long),'F140M native short and long exposures')
    result.sort(key=lambda r:(r['utc'],r['object']))
    with (ROOT/'perdita_history.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(result[0]));wr.writeheader();wr.writerows(result)
    fig,ax=plt.subplots(figsize=(9,3.6),layout='constrained')
    ax.axhline(0,color='#89929d',lw=.9,ls='--',zorder=1)
    for name,color,marker in [('Perdita','#155868','o'),('Belinda','#bd7351','x')]:
        data=[r for r in result if r['object']==name]
        xs=[date.fromisoformat(r['utc']) for r in data];ys=[r['longitude_residual_deg'] for r in data]
        ax.scatter(xs,ys,s=52,c=color,marker=marker,label=name,zorder=3)
        if name=='Perdita':
            for x,y in zip(xs,ys):ax.annotate(f'{y:.1f}°',(x,y),xytext=(0,-14 if y>-2 else 10),textcoords='offset points',ha='center',color=color,fontsize=9)
    ax.set_ylim(-23,4);ax.set_ylabel('Measured minus predicted\norbital longitude (degrees)')
    ax.xaxis.set_major_locator(dates.MonthLocator(bymonth=[1,7]));ax.xaxis.set_major_formatter(dates.DateFormatter('%b %Y'))
    ax.set_title('Perdita’s position relative to the JPL URA184 prediction',loc='left',fontweight='bold')
    ax.legend(frameon=False,loc='lower left',ncol=2);ax.grid(axis='y',alpha=.22);ax.spines[['top','right']].set_visible(False)
    fig.savefig(ROOT/'perdita_history.png',dpi=165)
if __name__=='__main__':main()
