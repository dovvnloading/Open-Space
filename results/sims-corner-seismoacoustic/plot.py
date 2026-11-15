"""Render waveforms and rank-preserving normalization from the verified cache."""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime

from analyze import ROOT, CENTERS, curves, trace_path, gate


def main():
    cmd = argparse.ArgumentParser(description=__doc__)
    cmd.add_argument('--cache',type=Path,required=True)
    cmd.add_argument('--output',type=Path,required=True)
    opt = cmd.parse_args()
    m=json.loads((ROOT/'manifest.json').read_text())
    fig,(left,right)=plt.subplots(1,2,figsize=(10.5,4.6),gridspec_kw={'width_ratios':[1.28,1]})
    colors={'ground':'#174a78','pressure':'#c36429'}
    stations=['SAW','EPH2','DY2','CHELN','CBS','WAT2','MNSN','NTAT']
    for i,sta in enumerate(stations):
        r=next(r for r in m['waveforms'] if r['event']=='wa' and r['nslc'].split('.')[1]==sta)
        a,_=curves(trace_path(opt.cache,r),UTCDateTime(m['sources']['wa']['cneos']['date']),(.6,6))
        mask=(CENTERS>=50)&(CENTERS<=275)
        norm=np.nanmax(a[mask])
        left.plot(CENTERS[mask],i+.66*a[mask]/norm,color=colors['ground'],lw=1.0)
    left.set(yticks=range(len(stations)),yticklabels=stations,xlabel='Seconds after optical peak',xlim=(50,275),
             title='Washington · near-field ground motion')
    alaska=['H21K','H22K','I21K','COLD']
    for i,sta in enumerate(alaska):
        for kind,chan in [('pressure','BDF'),('ground','BHZ')]:
            r=next(r for r in m['waveforms'] if r['event']=='ak' and r['nslc'].split('.')[1]==sta and r['nslc'].endswith(chan))
            a,_=curves(trace_path(opt.cache,r),UTCDateTime(m['sources']['ak']['cneos']['date']),(.6,6))
            if kind=='pressure':
                lo,hi=gate(r,35.5)
                peak=CENTERS[np.nanargmax(np.where((CENTERS>=lo)&(CENTERS<=hi),a,np.nan))]
            mask=(CENTERS>=peak-15)&(CENTERS<=peak+15)
            right.plot(CENTERS[mask]-peak,i+.75*a[mask]/np.nanmax(a[mask]),
                       color=colors[kind],lw=1.4,alpha=.86,label=kind.capitalize() if i==0 else None)
    right.set(yticks=range(4),yticklabels=alaska,xlabel='Seconds relative to local pressure peak',
              xlim=(-15,15),title='Alaska · colocated sensors')
    right.legend(loc='upper right',frameon=False,fontsize=9)
    for ax in (left,right):
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='x',alpha=.2)
    fig.tight_layout(pad=1.2)
    fig.savefig(opt.output,format='svg',metadata={'Date':None})


if __name__=='__main__':
    main()
