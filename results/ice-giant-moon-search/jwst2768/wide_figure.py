"""Summarize conditional F182M mosaic completeness and covered trial phases."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def wilson(k, n):
    p=k/n
    den=1+1.96**2/n
    center=(p+1.96**2/(2*n))/den
    half=1.96*np.sqrt(p*(1-p)/n+1.96**2/(4*n*n))/den
    return max(0,center-half),min(1,center+half)


def run(coverage, plants, destination):
    with Path(coverage).open() as f:
        grids=list(csv.DictReader(f))
    with Path(plants).open() as f:
        all_plants=list(csv.DictReader(f))
    kept=[r for r in all_plants if r['eligible']=='True']
    fields=[k for k in kept[0] if k.startswith('rank_')]
    rows=[]
    for key in fields:
        flux=float(key.removeprefix('rank_').removesuffix('uJy'))
        k=sum(float(r[key])>=1 for r in kept)
        low,high=wilson(k,len(kept))
        rows.append([flux,k,len(kept),k/len(kept),low,high])
    dest=Path(destination)
    with (dest/'completeness.csv').open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['uJy','recovered','eligible','fraction','wilson_low95','wilson_high95'])
        w.writerows(rows)
    with (dest/'completeness_by_radius.csv').open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['a_min_km','a_max_km','subset','draws','eligible','uJy','recovered','fraction','wilson_low95','wilson_high95'])
        regions=[(a,a+50000,'stratified') for a in range(200000,650000,50000)]
        regions.append((584000,650000,'beyond_Oberon'))
        for lo,hi,label in regions:
            group=[r for r in all_plants if lo<=float(r['a_km'])<hi]
            supported=[r for r in group if r['eligible']=='True']
            for key in fields:
                flux=float(key.removeprefix('rank_').removesuffix('uJy'))
                k=sum(float(r[key])>=1 for r in supported)
                n=len(supported)
                interval=wilson(k,n) if n else ('','')
                w.writerow([lo,hi,label,len(group),n,flux,k,k/n if n else '',*interval])
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,(ax,bx)=plt.subplots(1,2,figsize=(9.4,3.55),layout='constrained')
    radii=np.array([float(r['a_km']) for r in grids])/1000
    for key,name,color in [('supported_2visits','At least two visits','#136e7c'),('supported_3visits','All three visits','#a05237')]:
        frac=[int(r[key])/int(r['trials']) for r in grids]
        ax.plot(radii,frac,lw=1.8,label=name,color=color)
    ax.axvline(583.5,color='.5',ls=':',lw=1)
    ax.text(579,.58,'Oberon',rotation=90,ha='right',color='.3',fontsize=9)
    ax.set(xlabel='Circular orbit radius (thousand km)',ylabel='Fraction of trial phases with coverage',
           title='Three Uranus visits, F182M',ylim=(.45,1.01),xlim=(200,650))
    ax.legend(loc='lower left',frameon=False,fontsize=9)
    x,k,n,y,lo,hi=np.array(rows).T
    bx.fill_between(x,lo,hi,alpha=.2,color='#136e7c')
    bx.plot(x,y,'o-',color='#136e7c',lw=1.8,markersize=4)
    bx.axhline(.5,color='.5',ls=':',lw=1)
    bx.set(xlabel='Injected flux density (µJy)',ylabel='Fraction of supported plants recovered',
           title=f'Full-dither subset, {len(kept)} of {len(all_plants)} orbits',ylim=(-.02,1.03),xlim=(-.015,x.max()+.015))
    for a in (ax,bx):
        a.grid(axis='y',alpha=.18)
    fig.savefig(dest/'completeness.png',dpi=170)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('coverage');p.add_argument('plants');p.add_argument('destination')
    a=p.parse_args();run(a.coverage,a.plants,a.destination)
