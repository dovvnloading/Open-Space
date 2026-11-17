"""Plot signed 3.10–4.97 µm DR1 photometry with a fixed 0.12 µm display bin."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from frame_inputs import read_frames

import numpy as np

LABELS = {
    'WISE0855':'WISE 0855−0714 · known cold benchmark',
    'J041358':'J041358−475038 · published T9',
    'J181125':'J181125+665806 · published early-Y estimate',
    'J182532':'J182532+671906 · CatWISE trial track',
    'J085014':'J085014+541123 · CatWISE trial track',
    'J092244':'J092244+640537 · CatWISE trial track',
    'J010638.19+122908.3':'J010638.19+122908.3',
    'J213348.06+391040.9':'J213348.06+391040.9',
    'J192049.57+111222.6':'J192049.57+111222.6',
}
def main(inputs,output):
    entries=read_frames(inputs)
    panels=[(k,v) for k,v in LABELS.items() if any(r['target']==k for r in entries)]
    ncol=2 if len(panels)>3 else 1
    nrow=(len(panels)+ncol-1)//ncol
    fig,axes=plt.subplots(nrow,ncol,figsize=(5*ncol,2.8*nrow),sharex=True,squeeze=False,layout='constrained')
    for ax,(name,title) in zip(axes.flat,panels):
        rows=[r for r in entries if r['target']==name and r['status']=='fit' and 3.1<=float(r['wavelength_um'])<=4.97]
        if not rows:
            ax.set_visible(False)
            continue
        wave,flux,err=np.array([[float(r[k]) for k in ('wavelength_um','flux_ujy','sigma_ujy')] for r in rows]).T
        ax.axhline(0,color='0.55',lw=.8)
        for lower in np.arange(3.10,4.97,.12):
            m=(wave>=lower)&(wave<lower+.12)
            if not m.any():
                continue
            weight=err[m]**-2
            ax.errorbar(np.average(wave[m],weights=weight),np.average(flux[m],weights=weight),
                       yerr=1/np.sqrt(weight.sum()),fmt='o',ms=3.5,capsize=2,lw=1,color='#15588c')
        ax.set_title(title,loc='left',fontsize=10)
        ax.grid(axis='y',alpha=.2)
        ax.set_xlim(3.08,5.01)
    for ax in axes[:,0]:
        ax.set_ylabel('Flux density (µJy)')
    for ax in axes[-1]:
        ax.set_xlabel('Wavelength (µm)')
    title='SPHEREx DR1 trial-track photometry' if ncol==1 else 'SPHEREx DR1 flux measured at the predicted positions'
    fig.suptitle(title,fontsize=13 if ncol==1 else 15)
    fig.savefig(output,dpi=175)
if __name__=='__main__':
    main(sys.argv[1:-1],sys.argv[-1])
