"""SWAN diagnostic sensitivity study; not calibrated spectroscopy or discovery test."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
from astropy.io import fits
sys.path.insert(0,str(Path(__file__).resolve().parent))
from spherex_moving_screen import cache_path

def measure(im,var,flags,x,y,r,inner,outer):
 yy,xx=np.indices(im.shape);d=np.hypot(xx-x,yy-y)
 good=np.isfinite(im)&np.isfinite(var)&(var>0)&((flags&~(2097152|524288))==0)
 ap=d<=r;sky=(d>=inner)&(d<=outer)&good
 if not np.all(good[ap]) or sky.sum()<50:return None
 bg=np.median(im[sky]);err=np.sqrt(var[ap].sum()+ap.sum()**2*(1.2533*1.4826*np.median(abs(im[sky]-bg))/np.sqrt(sky.sum()))**2)
 return {'sum':float((im[ap]-bg).sum()),'formal_error':float(err),'sky_pixels':int(sky.sum())}

def main():
 rows=json.loads(Path(sys.argv[1]).read_text());out=[]
 for z in rows:
  if 'sha256' not in z:continue
  with fits.open(cache_path(Path('/tmp/spherex-cache'),z['url'])) as h:
   im=h['IMAGE'].data;v=h['VARIANCE'].data;f=h['FLAGS'].data;x=z['x'];y=z['y'];yy,xx=np.indices(im.shape)
   good=np.isfinite(im)&((f&~(2097152|524288))==0)
   rr=np.hypot(xx-x,yy-y);bg=np.median(im[(rr>=8)&(rr<=12)&good])
   # Local peak and positive-flux centroid are descriptive, not a PSF fit.
   roi=(rr<3)&good;w=np.where(roi,np.maximum(im-bg,0),0);cy=float((w*yy).sum()/w.sum());cx=float((w*xx).sum()/w.sum())
   i,j=np.unravel_index(np.argmax(np.where(roi,im,-np.inf)),im.shape)
   result={'observationid':z['observationid'],'start':z['start'],'approx_wavelength_um':z['wavelength_um'],'peak_offset_pixels':[float(j-x),float(i-y)],'positive_flux_centroid_offset_pixels':[cx-x,cy-y],'image_midpoint_mjd':float(h['IMAGE'].header['MJD-AVG']),'annulus_aperture_sensitivity':{f'{r}:{a}-{b}':measure(im,v,f,x,y,r,a,b) for r in [2.,2.5,3.5,5.] for a,b in [(6,10),(8,12),(12,17)]}}
   controls=[]
   for dx in [-18,-12,-6,0,6,12,18]:
    for dy in [-18,-12,-6,0,6,12,18]:
     if np.hypot(dx,dy)<10:continue
     if x+dx-2.5<0 or x+dx+2.5>=im.shape[1] or y+dy-2.5<0 or y+dy+2.5>=im.shape[0]:continue
     m=measure(im,v,f,x+dx,y+dy,2.5,6,10)
     if m:controls.append({'dx':dx,'dy':dy,**m})
   result['local_grid_controls']=controls;out.append(result)
 Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
