"""Compact-source aperture sensitivity, not a full comet coma correction."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates,gaussian_filter
from astropy.io import fits
from astropy.wcs import WCS
from astropy import log
log.setLevel('ERROR')
p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--cache',required=True);p.add_argument('--output',required=True);a=p.parse_args();out=[]
for row in json.load(open(a.input)):
 if not(row['start'].startswith('2026-01-02 10:') or row['start'].startswith('2025-12-25')):continue
 with fits.open(Path(a.cache)/(hashlib.sha256(row['url'].encode()).hexdigest()+'.fits')) as h:
  x,y=row['x'],row['y'];w=WCS(h['IMAGE'].header,key='A');w.sip=None;ax,ay=w.wcs_pix2world(x,y,0);tb=h['EPSF'].data;k=np.argmin((tb['XCENTER']-ax)**2+(tb['YCENTER']-ay)**2);ps=tb['EPSF'][k];yy,xx=np.indices(h['IMAGE'].data.shape);r=np.hypot(xx-x,yy-y);fractions={}
  for sigma in [0.,.5,1.]:
   # Zero padding avoids convolution-edge losses but cannot restore unmeasured wings.
   q=gaussian_filter(np.pad(ps,30),sigma*5) if sigma else np.pad(ps,30);q=q/q.sum();model=25*map_coordinates(q,[(yy-y)*5+46,(xx-x)*5+46],order=3,mode='constant',cval=0);b=np.median(model[(r>=8)&(r<=12)]);fractions[str(sigma)]={str(rad):float((model[r<=rad]-b).sum()) for rad in [2.5,3.5,5.]}
  out.append(dict(observationid=row['observationid'],start=row['start'],wavelength_approx=row['wavelength_um'],epsf_index=int(k),fractions=fractions,measured=row['diagnostic_outlier_inclusive']))
Path(a.output).write_text(json.dumps(out,indent=2)+'\n')
for x in out:print(x['start'],x['wavelength_approx'],[x['fractions'][str(s)]['2.5'] for s in [0.,.5,1.]])
