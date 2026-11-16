"""Inspect every LoDeSS catalogue component within 45 arcsec of ILTJ100005.91+812702.7.

python J100005_images.py DATA_DIR OUTPUT_DIR
Requires run of map_fields.py. Acquires the required images from the official releases.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from matplotlib.patches import Circle

from acquire_lodess import retrieve, digest
from photometry import RadioImage
from resolution_check import aperture

def run(root,out):
    out.mkdir(parents=True,exist_ok=True)
    sky=SkyCoord(150.02463*u.deg,81.45076*u.deg)
    fields=pd.read_csv(root/'fields.csv').set_index('field')
    rows=[]
    for tab, mapping in [('sources','source_fields'),('gaussians','gaussian_fields')]:
        data=fits.getdata(root/(tab+'.fits'))
        names=pd.read_csv(root/(mapping+'.csv')).field
        distance=sky.separation(SkyCoord(data['RA']*u.deg,data['DEC']*u.deg)).arcsec
        for i in np.flatnonzero(distance<45):
            row=dict(catalogue=tab,row=i,field=names.iloc[i],sep_as=distance[i])
            for c in ['Source_id','Isl_id','Gaus_id','RA','DEC','Total_flux','E_Total_flux','Peak_flux','Isl_rms','S_Code']:
                if c in data.names:row[c]=data[c][i]
            rows.append(row)
    pd.DataFrame(rows).to_csv(out/'components.csv',index=False)
    lodessfields=sorted({r['field'] for r in rows})
    items=[]
    for field in lodessfields:
        url=f'https://lofar-surveys.org/public/LoDeSS/{field}/{field}-MFS-I-image-pb.fits'
        items.append((field, root/'images'/(field+'-I.fits'),url))
    for beam in [6,20]:
        suffix='-low' if beam==20 else ''
        pos=sky.to_string(style='hmsdms',sep=':',precision=4)
        url=f'https://lofar-surveys.org/dr3{suffix}-cutout.fits?'+urlencode({'pos':pos,'size':12})
        items.append((f'LoTSS-{beam}',out/f'LoTSS-{beam}.fits',url))
    def get(item):
        name,f,url=item
        if not f.exists():retrieve(url,f)
        return dict(product=name,sha256=digest(f),url=url)
    with ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(get,items))
    pd.DataFrame(records).to_csv(out/'manifest.csv',index=False)
    fig,axes=plt.subplots(1,len(items),figsize=(2.95*len(items),3.7))
    measurements=[]
    for ax,(name,f,url) in zip(axes,items):
        with RadioImage(f) as im:
            h=im.hdus[0].header
            forced=im.measure(sky.ra.deg,sky.dec.deg) or {}
            cut=Cutout2D(im.data,sky,(7*u.arcmin,7*u.arcmin),wcs=im.wcs)
            x,y=cut.wcs.world_to_pixel(sky)
            beam=3600*h['BMAJ']
            flux,bg,noise=aperture(cut.data,x,y,im.pixel_arcsec,beam)
            row=dict(image=name,date_obs=h.get('DATE-OBS'),beam_arcsec=beam,
                     ap60_jy=flux,pixel_rms_jy_beam=noise,**forced)
            if name in fields.index:
                q=fields.loc[name]
                row.update(field_quality=q.quality,
                           offaxis_deg=sky.separation(SkyCoord(q.ra_deg*u.deg,q.dec_deg*u.deg)).deg)
            measurements.append(row)
            hdr=cut.wcs.to_header()
            for key in ['BMAJ','BMIN','BPA','BUNIT','DATE-OBS']:
                if key in h:hdr[key]=h[key]
            fits.writeto(out/(name+'-stamp.fits'),cut.data,hdr,overwrite=True)
            stamp=Cutout2D(cut.data-bg,sky,(2.5*u.arcmin,2.5*u.arcmin),wcs=cut.wcs)
            x,y=stamp.wcs.world_to_pixel(sky)
            ax.imshow(stamp.data,origin='lower',cmap='inferno',interpolation='nearest',
                      vmin=-2*noise,vmax=max(4*noise,np.nanmax(stamp.data)*.7))
            ax.add_patch(Circle((x,y),22.5/im.pixel_arcsec,fill=False,color='cyan',lw=.8))
            title=name if name in fields.index else 'LoTSS 144 MHz'
            line=(h.get('DATE-OBS','')[:10]+' · 23 MHz') if name in fields.index else f'{beam:.0f}″ beam'
            ax.set_title(title+'\n'+line,fontsize=9)
            ax.set_xticks([]);ax.set_yticks([])
    pd.DataFrame(measurements).to_csv(out/'image_photometry.csv',index=False)
    fig.suptitle('ILTJ100005.91+812702.7 · 2.5 arcmin field · separate brightness scales\nCyan circles show the 45″ LoDeSS beam at the LoTSS position')
    fig.tight_layout()
    fig.savefig(out/'cutouts.png',dpi=150)
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();run(a.root,a.out)
