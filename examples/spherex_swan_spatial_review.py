"""Independent spatial diagnostic using same FITS data. Not independent observation."""
import json,hashlib,argparse
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates,gaussian_filter
from scipy.optimize import least_squares
parser=argparse.ArgumentParser();parser.add_argument('input');parser.add_argument('output');parser.add_argument('--cache',default='/tmp/spherex-cache');args=parser.parse_args()
rows=json.load(open(args.input));results=[]
for row in rows:
 p=Path(args.cache)/(hashlib.sha256(row['url'].encode()).hexdigest()+'.fits')
 with fits.open(p) as h:
  im=h['IMAGE'].data;var=h['VARIANCE'].data;flags=h['FLAGS'].data
  x,y=WCS(h['IMAGE'].header,h).all_world2pix([[row['ra'],row['dec']]],0)[0]
  ax,ay=WCS(h['IMAGE'].header,h,key='A').wcs_pix2world([[x,y]],0)[0]
  tb=h['EPSF'].data;k=int(np.argmin((tb['XCENTER']-ax)**2+(tb['YCENTER']-ay)**2));ps=tb['EPSF'][k]; yy,xx=np.indices(im.shape)
  good=np.isfinite(im)&np.isfinite(var)&(var>0)&((flags&~(2097152|524288))==0);r=np.hypot(xx-x,yy-y)
  out={'observationid':row['observationid'],'start':row['start'],'url':row['url'],'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'wavelength_approx':row.get('wavelength_um'),'epsf_index':k,'epsf_center':[float(tb['XCENTER'][k]),float(tb['YCENTER'][k])],'detector_position':[float(ax),float(ay)],'target_pixel':[float(x),float(y)],'invalid_core_pixels':int(((r<=2.5)&~good).sum())}
  if out['invalid_core_pixels']:
   out['status']='rejected core invalid flags';results.append(out);continue
  # radius2.5: fits within33/5 native-pixel PSF support at small shifts.
  mask=(r<=2.5)&good;X=xx[mask];Y=yy[mask];z=im[mask];s=np.sqrt(var[mask]);sky=(r>8)&(r<12)&good;bg=float(np.median(im[sky]))
  def model(v):
   amp,b,dx,dy=v;return amp*map_coordinates(ps,[(Y-y-dy)*5+16,(X-x-dx)*5+16],order=3,mode='constant',cval=0)*25+b
  fit=least_squares(lambda v:(model(v)-z)/s,[float(np.sum(z-bg)),bg,0,0],bounds=([0,-np.inf,-.5,-.5],[np.inf,np.inf,.5,.5]))
  out['fit']={'amp':float(fit.x[0]),'background':float(fit.x[1]),'offset':[float(v) for v in fit.x[2:]],'chi2':float(np.sum(fit.fun**2)),'dof':int(len(z)-4),'reduced_chi2':float(np.sum(fit.fun**2)/(len(z)-4)),'n':len(z)}
  # Nonparametric growth, explicit background and positivity centroid.
  out['growth']={str(rad):float(np.sum(im[(r<=rad)&good]-bg)) for rad in [1,1.5,2,2.5,3,4,5]}
  q=np.where((r<=3)&good,np.maximum(im-bg,0),0);out['positive_centroid_offset']=[float((xx*q).sum()/q.sum()-x),float((yy*q).sum()/q.sum()-y)]
  out['status']='diagnostic OUTLIER-inclusive, not calibrated flux';results.append(out)
Path(args.output).write_text(json.dumps(results,indent=2))
for r in results:print(r['start'],r['wavelength_approx'],r.get('fit'),r.get('positive_centroid_offset'),r.get('status'))
