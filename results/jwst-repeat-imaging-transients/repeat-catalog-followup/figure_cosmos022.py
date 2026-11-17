"""Render the image evidence and constituent-exposure light curve of COSMOS-022."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from campaign import HERE


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stamp_cache',type=Path)
    cache=p.parse_args().stamp_cache
    key='cosmos_medium-022'
    pair=pd.read_csv(HERE/'pairs.csv').query("pair=='cosmos_medium' and band=='f162m'").iloc[0]
    scale=pd.read_csv(HERE/'screen_counts.csv').query("pair=='cosmos_medium' and band=='f162m'").flux_scale.iloc[0]
    a,b=[np.load(cache/f'{key}_{i}.npz')['pixels'][25:78,25:78] for i in (pair.reference,pair.search)]
    a=a*scale;difference=b-a
    fig=plt.figure(figsize=(9,5.9),layout='constrained')
    grid=fig.add_gridspec(2,3,height_ratios=(1,1.2))
    fig.suptitle('COSMOS-022   |   NIRCam repeat imaging',fontsize=15,fontweight='bold')
    bounds=np.nanpercentile(np.concatenate([a.ravel(),b.ravel()]),[4,99.7])
    norm=AsinhNorm(linear_width=(bounds[1]-bounds[0])*.09,vmin=bounds[0],vmax=bounds[1])
    residual=np.nanpercentile(abs(difference),99)
    for j,(im,label) in enumerate([(a,'F162M  ·  Dec–Jan'),(b,'F162M  ·  May'),(difference,'May minus Dec–Jan')]):
        ax=fig.add_subplot(grid[0,j])
        if j==2:ax.imshow(im,origin='lower',cmap='RdBu_r',vmin=-residual,vmax=residual,interpolation='nearest')
        else:ax.imshow(im,origin='lower',cmap='gray',norm=norm,interpolation='nearest')
        ax.set(xticks=[],yticks=[],title=label)
        ax.plot([4,16.5],[4,4],lw=2,color='#e0b94e')
        ax.text(4,6,'0.5″',color='#f0c559',fontsize=8)
    ax=fig.add_subplot(grid[1,:]);t=pd.read_csv(HERE/'individual_apertures.csv').query('candidate==@key')
    dates={61038:'Dec 29, 2025',61044:'Jan 4, 2026',61178:'May 18, 2026',61184:'May 24, 2026'}
    groups=sorted(t.mjd.astype(int).unique())
    for k,(band,color) in enumerate([('f162m','#156963'),('f182m','#a54a21')]):
        for j,day in enumerate(groups):
            dat=t[(t.band==band)&(t.mjd.astype(int)==day)]
            x=j+(k-.5)*.23
            ax.scatter(x+np.linspace(-.045,.045,len(dat)),dat.fixed_njy,color=color,s=18,alpha=.8,label=band.upper() if j==0 else None)
            ax.hlines(np.median(dat.fixed_njy),x-.082,x+.082,color=color,lw=2.5)
    ax.set_xticks(range(len(groups)),[dates.get(int(d),f'MJD {int(d)}') for d in groups])
    ax.set(ylabel='Flux in 0.16″ aperture (nJy)',title='Individual exposures at the same sky position')
    ax.legend(ncols=2,frameon=False,loc='upper left')
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.22)
    fig.savefig(HERE/'cosmos022.png',dpi=150)
    plt.close(fig)


if __name__=='__main__':main()
