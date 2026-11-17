"""Larger region connecting a multi-dither NISP/VIS linear feature to a bright Gaia star."""
import argparse
import csv
import json
from io import BytesIO
from pathlib import Path
from urllib.parse import urlencode

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

from acquire import IBE,fetch,sha

ID='102158277_positive_35'
STAR=(273.80322781496153,64.49441035233595)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('candidates','products','cache','out'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();r=next(z for z in csv.DictReader(a.candidates.open()) if z['id']==ID)
    product=next(z for z in csv.DictReader(a.products.open()) if z['band_or_product']=='VIS' and z['tile']==r['tile'])
    mid=np.array([float(r['ra_deg']),float(r['dec_deg'])]);ctr=(mid+STAR)/2
    url=IBE+'q1/MER/'+product['relative_to_mer']+f'?center={ctr[0]:.8f},{ctr[1]:.8f}deg&size=140arcsec&gzip=false'
    raw=fetch(url,a.cache/'context_102158277.fits')
    h=fits.getheader(BytesIO(raw));im=fits.getdata(BytesIO(raw));w=WCS(h)
    star,center=w.all_world2pix([STAR,mid],0);ends=w.all_world2pix([[float(r[f'ra{i}_deg']),float(r[f'dec{i}_deg'])] for i in (0,1)],0)
    direction=(center-star)/np.linalg.norm(center-star);normal=np.array([-direction[1],direction[0]])
    sx,sy=np.meshgrid(np.linspace(-16,118,1340),np.linspace(-18,18,360))
    pix=star+sx[:,:,None]*direction/.1+sy[:,:,None]*normal/.1
    img=map_coordinates(im,[pix[:,:,1],pix[:,:,0]],order=1,mode='constant',cval=np.nan)
    med=np.nanmedian(img);sig=np.nanmedian(np.abs(img-med))*1.4826
    f,ax=plt.subplots(figsize=(10,3.25))
    ax.imshow(np.arcsinh((img-med)/(2*sig)),cmap='gray',vmin=-1,vmax=3.2,extent=[-16,118,-18,18],aspect='equal',origin='lower')
    endpos=np.sort((ends-star)@direction*.1)
    ax.plot(endpos,[4,4],color='#53d8fc',lw=2);ax.plot(endpos,[4,4],linestyle='none',marker='|',color='#53d8fc',ms=7)
    ax.annotate('Gaia G = 10.52',xy=(0,0),xytext=(5,12),color='white',arrowprops=dict(arrowstyle='-',color='white'))
    ax.text(np.mean(endpos),6,'18″ detector interval',ha='center',color='#53d8fc')
    ax.set_xlabel('Offset along the star–feature direction (arcsec)');ax.set_ylabel('Cross-axis offset (arcsec)')
    f.tight_layout();f.savefig(a.out/'stellar_spike_context.jpg',dpi=150);plt.close(f)
    (a.out/'stellar_spike_context.json').write_text(json.dumps(dict(id=ID,star_id=r['gaia_dr3_spike_id'],star_ra_deg=STAR[0],star_dec_deg=STAR[1],star_g=r['gaia_dr3_g'],
        data_url=url,bytes=len(raw),sha256=sha(raw),star_to_line_center_arcsec=np.linalg.norm(center-star)*.1,
        star_cross_track_arcsec=r['gaia_dr3_cross_arcsec']),indent=1)+'\n')


if __name__=='__main__':main()
