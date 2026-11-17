"""Summarize complete-exposure-block dropout and morphology of the five Uranus cells."""
import csv
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from coverage import OUT
from miri_audit import write_csv

if __name__=='__main__':
    rows=list(csv.DictReader((OUT/'uranus_native_vet.csv').open()));summary=[]
    stamps=list(csv.DictReader((OUT/'uranus_stamps.csv').open()));shapes=[]
    bands=['F560W','F1800W','F1000W'];fig,axes=plt.subplots(1,3,figsize=(8.3,3.1),sharex=True,sharey=True,layout='constrained')
    for j,v in enumerate(sorted({r['visit'] for r in rows})):
        for template in sorted({r['template'] for r in rows}):
            use=[r for r in rows if r['usable']=='1' and r['template']==template and r['visit']==v]
            for block in sorted({r['file'].split('_')[1] for r in use}):
                z=np.clip([float(r['normalized_score']) for r in use if r['file'].split('_')[1]!=block],-7,7)
                summary.append(dict(template=template,visit=v,filter=bands[j],dropped_activity=block,remaining=len(z),normalized_remaining=z.sum()/np.sqrt(len(z))))
        array=np.array([float(r['median_normalized_score']) for r in stamps if r['visit']==v]).reshape(25,25)
        yy,xx=np.mgrid[-12:13,-12:13];loc=np.unravel_index(np.nanargmax(np.where(np.hypot(xx,yy)<=10,array,np.nan)),array.shape)
        # Remote on-row samples measure whether the response extends far outside one PSF across the detector.
        wing=(np.abs(yy)<=1)&(np.abs(xx)>=8)&(np.abs(xx)<=12)
        shapes.append(dict(visit=v,filter=bands[j],score_at_prediction=array[12,12],peak_in_10px=array[loc],peak_dx=xx[loc],peak_dy=yy[loc],
                           median_remote_same_row=np.nanmedian(array[wing]),max_remote_same_row=np.nanmax(array[wing])))
        im=axes[j].imshow(array,origin='lower',extent=(-12.5,12.5,-12.5,12.5),cmap='viridis',vmin=-2.5,vmax=2.5)
        axes[j].plot(0,0,'+',color='white',ms=10);axes[j].set_title(bands[j]);axes[j].set_xlabel('Detector offset (pixels)')
    axes[0].set_ylabel('Detector offset (pixels)');fig.colorbar(im,ax=axes,shrink=.65,label='Median native score')
    fig.savefig(OUT/'uranus_source_stamps.png',dpi=150);plt.close(fig)
    write_csv(OUT/'uranus_activity_dropout.csv',list(summary[0]),summary);write_csv(OUT/'uranus_morphology.csv',list(shapes[0]),shapes)
