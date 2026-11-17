"""Examine wide-field context around the projected galaxy/line intersection."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from PIL import Image

from targeted import fetch
HERE=Path(__file__).resolve().parent
OID='2699220004635136134'

def run():
    keys={(r['tile'],r['band_or_product']):r['relative_to_mer'] for r in csv.DictReader((HERE.parent/'euclid/q1_products.csv').open())}
    r=next(r for r in csv.DictReader((HERE/'hats_selected.csv').open()) if r['object_id']==OID)
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,5),layout='constrained')
    records=[]
    for band,ax in zip(['VIS','H'],axes):
        img,x,y,rec=fetch(keys[(r['tileid'],band)],r,180)
        ok=np.isfinite(img)&(img!=0)
        med,q=np.percentile(img[ok],[50,84])
        view=(img-med)/max(q-med,1e-10)
        ax.imshow(np.where(ok,np.arcsinh(view),np.nan),vmin=-1.5,vmax=4.7,origin='lower',cmap='gray',
                  extent=[(-x-.5)*.1,(img.shape[1]-x-.5)*.1,(-y-.5)*.1,(img.shape[0]-y-.5)*.1],interpolation='nearest')
        ax.plot([0,0],[-11,-6],color='darkorange',lw=1)
        ax.plot([-11,-6],[0,0],color='darkorange',lw=1)
        ax.set(xlabel='Mosaic x offset (arcsec)',title=band)
        records.append(dict(band=band,nx=img.shape[1],ny=img.shape[0],valid_fraction=ok.mean(),**rec))
    axes[0].set_ylabel('Mosaic y offset (arcsec)')
    fig.savefig(HERE/'galaxy_projection.jpg',dpi=155)
    (HERE/'galaxy_projection.json').write_text(json.dumps({'id':OID,'ra':r['ra'],'dec':r['dec'],'images':records},indent=2)+'\n')
if __name__=='__main__':run()
