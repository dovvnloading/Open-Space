"""Test the catalogue line that cuts across a galaxy against a distant Gaia star."""
import csv
import io
import json
import hashlib
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates

HERE=Path(__file__).resolve().parent
OID='2699220004635136134'
STAR='2160979470332604672'
def run():
    r=next(r for r in csv.DictReader((HERE/'hats_selected.csv').open()) if r['object_id']==OID)
    k=next(k['relative_to_mer'] for k in csv.DictReader((HERE.parent/'euclid/q1_products.csv').open()) if k['tile']==r['tileid'] and k['band_or_product']=='VIS')
    point=SkyCoord(float(r['ra'])*u.deg,float(r['dec'])*u.deg)
    star=SkyCoord(270.06650792245415*u.deg,63.55150353883631*u.deg)
    pa=point.position_angle(star);sep=point.separation(star)
    cen=point.directional_offset_by(pa,sep/2)
    url='https://irsa.ipac.caltech.edu/ibe/data/euclid/q1/MER/'+k+f'?center={cen.ra.deg:.8f},{cen.dec.deg:.8f}deg&size=360arcsec&gzip=false'
    with urlopen(url,timeout=180) as f:b=f.read()
    with fits.open(io.BytesIO(b)) as h:
        im=h[0].data;w=WCS(h[0].header)
        x0,y0=w.world_to_pixel(point);x1,y1=w.world_to_pixel(star)
        xx=np.array([x0,x1]);yy=np.array([y0,y1])
        radial=np.hypot(x1-x0,y1-y0);unit=np.array([x1-x0,y1-y0])/radial
        fig,ax=plt.subplots(figsize=(7.5,7.2),layout='constrained')
        sm=im[::3,::3]
        mask=np.isfinite(sm)&(sm!=0);med,q=np.percentile(sm[mask],[50,84])
        ax.imshow(np.where(mask,np.arcsinh((sm-med)/max(q-med,1e-9)),np.nan),origin='lower',vmin=-1,vmax=6.2,cmap='gray',
                  extent=[(-x0-1.5)*.1,(sm.shape[1]*3-x0-1.5)*.1,(-y0-1.5)*.1,(sm.shape[0]*3-y0-1.5)*.1])
        for x,y,label,offset in [(0,0,'MER line',(-25,-25)),((x1-x0)*.1,(y1-y0)*.1,'Gaia G = 6.75',(25,25))]:
            ax.annotate(label,xy=(x,y),xytext=offset,textcoords='offset points',color='#ffa82d',ha='center',fontsize=10,
                        arrowprops=dict(arrowstyle='->',color='#ffa82d',shrinkB=4))
        ax.set(xlabel='Mosaic x offset (arcsec)',ylabel='Mosaic y offset (arcsec)',title='Euclid VIS · catalogue line and stellar diffraction pattern')
        fig.savefig(HERE/'spike_context.jpg',dpi=145)
        # Record the peak offset of stacked transverse profiles in ten equal axial intervals.
        v=np.array([-unit[1],unit[0]]);cross=np.arange(-40,41);parts=[]
        for i in range(10):
            t=np.linspace((.01+.098*i)*radial,(.01+.098*(i+1))*radial,101)
            x=x0+unit[0]*t+v[0]*cross[:,None];y=y0+unit[1]*t+v[1]*cross[:,None]
            traces=map_coordinates(im,[y,x],order=1,mode='constant',cval=np.nan)
            traces[traces==0]=np.nan
            z=np.nanmedian(traces,axis=1)
            if not np.isfinite(z).any():parts.append(None);continue
            base=float(np.nanmedian(z[abs(cross)>22]));peak=np.nanargmax(z)
            parts.append(dict(bin=i+1,offset_arcsec=float(cross[peak]*.1),contrast=float(z[peak]-base)))
    (HERE/'spike_context.json').write_text(json.dumps(dict(object_id=OID,gaia_source=STAR,gmag=6.7483587,star_icrs_deg=[star.ra.deg,star.dec.deg],
        separation_arcsec=sep.arcsec,star_pa_deg=pa.deg,mer_position_angle_deg=float(r['mer_position_angle']),
        url=url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),profiles=parts),indent=2)+'\n')
if __name__=='__main__':run()
