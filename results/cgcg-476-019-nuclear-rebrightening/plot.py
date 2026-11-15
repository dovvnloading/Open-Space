"""Draw the light curves used in the report. Requires matplotlib."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from verify import ROOT, TARGET, DISCOVERY_MJD, direct_photometry, infrared_visits

def year(mjd):
    return 2000 + (np.asarray(mjd)-51544)/365.2425

def main():
    visits=infrared_visits()
    target=visits[visits.source==TARGET]
    optical=direct_photometry()
    fig,axes=plt.subplots(2,1,figsize=(9.4,6.5),sharex=True,constrained_layout=True)
    for filt,label,color in [('zg','ZTF g','#087d79'),('zr','ZTF r','#9e3d20')]:
        x=optical[(optical.field==647)&(optical.filtercode==filt)].copy()
        x['bin']=(x.mjd//45).astype(int)
        bins=x.groupby('bin').agg(mjd=('mjd','median'),mag=('mag','median'),n=('mjd','size'))
        bins=bins[bins.n>=4]
        baseline=x[x.year<2021].mag.median()
        axes[0].scatter(year(bins.mjd),baseline-bins.mag,s=19,color=color,label=label)
    for band,label,c in [('w1','W1 (3.4 μm)','#8561a2'),('w2','W2 (4.6 μm)','#18538b')]:
        axes[1].plot(year(target.mjd),target[band+'_excess_mag'],'o-',color=c,lw=1.3,ms=4.8,label=label)
    for ax in axes:
        ax.axhline(0,c='0.7',lw=0.9)
        ax.axvline(year(DISCOVERY_MJD),c='0.35',ls='--',lw=1)
        ax.spines[['right','top']].set_visible(False)
        ax.legend(frameon=False,loc='upper left',ncols=2)
        ax.set_ylabel('Brightening (mag)')
        ax.grid(axis='y',color='0.92')
    axes[0].text(year(DISCOVERY_MJD)+.1,-.06,'SN 2021ued',size=9)
    axes[0].set_title('Optical brightness of the galaxy nucleus',loc='left')
    axes[1].set_title('Infrared brightness relative to the same observing season',loc='left')
    axes[1].set_xlabel('Calendar year'); axes[1].set_xticks(range(2014,2027,2))
    fig.savefig(ROOT/'light-curves.png',dpi=155,bbox_inches='tight')
    plt.close(fig)
if __name__=='__main__':
    main()
